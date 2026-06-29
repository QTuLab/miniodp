#!/usr/bin/env python3
"""
Build a miniodp geneinfo database from local annotation files.

Supported inputs:
- GTF/GFF-like annotation with gene_id attributes
- Optional gene-level annotation TSV with gene symbol and description columns
- Optional gene-to-GO/KEGG TSV file
- Optional GO ontology file for GO term names and namespaces
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import logging
import re
import sqlite3
import tempfile
import time
import urllib.request
import zipfile
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from logging_utils import configure_script_logger, create_print_proxy


logger = configure_script_logger(__name__)
print = create_print_proxy(logger)

GO_BASIC_OBO_URLS = [
    "https://current.geneontology.org/ontology/go-basic.obo",
    "https://purl.obolibrary.org/obo/go/go-basic.obo",
]

ATTR_PATTERN = re.compile(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"')
MISSING_TEXT_VALUES = {"", "na", "n/a", "-", ".", "null", "none"}


@contextmanager
def open_text_auto(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield handle
        return
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            members = [info.filename for info in archive.infolist() if not info.is_dir()]
            if len(members) != 1:
                raise ValueError(f"Expected exactly one file inside zip archive: {path}")
            with archive.open(members[0], "r") as raw_handle:
                with io.TextIOWrapper(raw_handle, encoding="utf-8") as text_handle:
                    yield text_handle
        return
    with open(path, "r", encoding="utf-8") as handle:
        yield handle


def parse_attributes(raw: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for key, value in ATTR_PATTERN.findall(raw):
        attrs[key] = value

    if attrs:
        return attrs

    for item in raw.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def extract_go_definition(raw_value: str) -> str:
    raw_value = raw_value.strip()
    if raw_value.startswith('"') and '"' in raw_value[1:]:
        return raw_value.split('"', 2)[1]
    return raw_value


def normalize_optional_text(raw_value: Optional[str]) -> str:
    value = (raw_value or "").strip()
    if value.lower() in MISSING_TEXT_VALUES:
        return ""
    return value


def normalize_header_name(raw_value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", raw_value.strip().lower()).strip("_")


class LocalGeneinfoBuilder:
    def __init__(
        self,
        *,
        gtf_file: Path,
        output_dir: Path,
        annotation_tsv: Optional[Path] = None,
        go_kegg_tsv: Optional[Path] = None,
        go_obo: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        download_go_obo: bool = False,
        species_label: str = "local_species",
    ) -> None:
        self.gtf_file = gtf_file
        self.annotation_tsv = annotation_tsv
        self.go_kegg_tsv = go_kegg_tsv
        self.go_obo = go_obo
        self.cache_dir = cache_dir
        self.download_go_obo = download_go_obo
        self.species_label = species_label
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.geneinfo_db = self.output_dir / "geneinfo.db"
        self.geneinfo_toml = self.output_dir / "geneinfo.toml"

    def validate_inputs(self) -> None:
        if not self.gtf_file.exists():
            raise FileNotFoundError(f"GTF file not found: {self.gtf_file}")
        if self.annotation_tsv and not self.annotation_tsv.exists():
            raise FileNotFoundError(f"Annotation TSV file not found: {self.annotation_tsv}")
        if self.go_kegg_tsv and not self.go_kegg_tsv.exists():
            raise FileNotFoundError(f"GO/KEGG TSV file not found: {self.go_kegg_tsv}")
        if self.go_obo and not self.go_obo.exists():
            raise FileNotFoundError(f"GO ontology file not found: {self.go_obo}")

    def load_annotation_map(self) -> Dict[str, Dict[str, str]]:
        if not self.annotation_tsv:
            print("Annotation TSV not provided; ENS_INFO names and descriptions will come from GTF only")
            return {}

        print(f"Parsing gene annotation from {self.annotation_tsv}")
        annotation_map: Dict[str, Dict[str, str]] = {}

        with open_text_auto(self.annotation_tsv) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError(f"Annotation TSV has no header: {self.annotation_tsv}")

            header_map = {
                normalize_header_name(field_name): field_name
                for field_name in reader.fieldnames
                if field_name is not None
            }

            gene_id_key = next(
                (
                    header_map[key]
                    for key in ("geneid", "gene_id", "ens", "gene")
                    if key in header_map
                ),
                None,
            )
            if not gene_id_key:
                raise ValueError(
                    "Annotation TSV must include a gene ID column such as geneID, gene_id, or ENS"
                )

            gene_name_key = next(
                (
                    header_map[key]
                    for key in ("gene_symbol", "genesymbol", "gene_name", "genename", "symbol")
                    if key in header_map
                ),
                None,
            )
            description_key = next(
                (
                    header_map[key]
                    for key in ("description", "desc", "annotation", "product")
                    if key in header_map
                ),
                None,
            )

            gene_name_count = 0
            description_count = 0

            for row in reader:
                gene_id = normalize_optional_text(row.get(gene_id_key))
                if not gene_id:
                    continue

                gene_name = normalize_optional_text(row.get(gene_name_key)) if gene_name_key else ""
                description = normalize_optional_text(row.get(description_key)) if description_key else ""

                if gene_name:
                    gene_name_count += 1
                if description:
                    description_count += 1

                previous = annotation_map.get(gene_id)
                if previous:
                    if gene_name and previous["gene_name"] and gene_name != previous["gene_name"]:
                        raise ValueError(f"Conflicting gene names for {gene_id}: {previous['gene_name']} vs {gene_name}")
                    if description and previous["description"] and description != previous["description"]:
                        raise ValueError(
                            f"Conflicting descriptions for {gene_id}: {previous['description']} vs {description}"
                        )
                    if gene_name and not previous["gene_name"]:
                        previous["gene_name"] = gene_name
                    if description and not previous["description"]:
                        previous["description"] = description
                    continue

                annotation_map[gene_id] = {
                    "gene_name": gene_name,
                    "description": description,
                }

        print(
            f"Prepared {len(annotation_map):,} annotation rows "
            f"({gene_name_count:,} with gene symbols, {description_count:,} with descriptions)"
        )
        return annotation_map

    def parse_gtf(self) -> Tuple[List[Tuple[str, str, str, str]], List[Tuple[str, str, int, int, int]]]:
        print(f"Parsing annotation from {self.gtf_file}")

        gene_records: Dict[str, Dict[str, object]] = {}

        with open_text_auto(self.gtf_file) as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) < 9:
                    continue

                chrom, _source, feature, start, end, _score, strand_symbol, _phase, attrs_raw = parts
                attrs = parse_attributes(attrs_raw)
                gene_id = (
                    attrs.get("gene_id")
                    or attrs.get("gene")
                    or attrs.get("ID")
                )
                if not gene_id:
                    continue

                gene_name = (
                    attrs.get("gene_name")
                    or attrs.get("Name")
                    or attrs.get("gene_symbol")
                    or gene_id
                )
                gene_type = (
                    attrs.get("gene_biotype")
                    or attrs.get("gene_type")
                    or attrs.get("transcript_biotype")
                    or attrs.get("transcript_type")
                    or ""
                )
                description = (
                    attrs.get("description")
                    or attrs.get("product")
                    or attrs.get("Note")
                    or attrs.get("note")
                    or ""
                )

                start_int = int(start)
                end_int = int(end)
                strand = 1 if strand_symbol == "+" else -1 if strand_symbol == "-" else 0

                record = gene_records.setdefault(
                    gene_id,
                    {
                        "gene_name": gene_name,
                        "gene_type": gene_type,
                        "description": description,
                        "chromosome": chrom,
                        "start": start_int,
                        "end": end_int,
                        "strand": strand,
                        "has_cds": False,
                    },
                )

                record["start"] = min(record["start"], start_int)
                record["end"] = max(record["end"], end_int)

                if not record["gene_name"] and gene_name:
                    record["gene_name"] = gene_name
                if not record["gene_type"] and gene_type:
                    record["gene_type"] = gene_type
                if not record["description"] and description:
                    record["description"] = description
                if feature == "CDS":
                    record["has_cds"] = True

                if record["chromosome"] != chrom:
                    raise ValueError(f"Gene {gene_id} spans multiple chromosomes: {record['chromosome']} vs {chrom}")
                if record["strand"] not in (0, strand):
                    raise ValueError(f"Gene {gene_id} spans multiple strands")

        annotation_map = self.load_annotation_map()
        missing_genes = sorted(gene_id for gene_id in annotation_map if gene_id not in gene_records)
        if missing_genes:
            preview = ", ".join(missing_genes[:10])
            raise ValueError(
                f"Annotation TSV contains {len(missing_genes)} gene IDs not present in GTF: {preview}"
            )

        updated_gene_names = 0
        updated_descriptions = 0
        for gene_id, annotation in annotation_map.items():
            record = gene_records[gene_id]

            gene_name = annotation["gene_name"]
            if gene_name and gene_name != record["gene_name"]:
                record["gene_name"] = gene_name
                updated_gene_names += 1

            description = annotation["description"]
            if description and description != record["description"]:
                record["description"] = description
                updated_descriptions += 1

        if annotation_map:
            print(
                f"Applied annotation TSV to ENS_INFO "
                f"({updated_gene_names:,} gene names updated, {updated_descriptions:,} descriptions updated)"
            )

        info_rows: List[Tuple[str, str, str, str]] = []
        locus_rows: List[Tuple[str, str, int, int, int]] = []

        for gene_id, record in sorted(gene_records.items()):
            gene_type = record["gene_type"] or ("protein_coding" if record["has_cds"] else "unknown")
            info_rows.append(
                (
                    gene_id,
                    str(record["gene_name"] or gene_id),
                    str(gene_type),
                    str(record["description"] or ""),
                )
            )
            locus_rows.append(
                (
                    gene_id,
                    str(record["chromosome"]),
                    int(record["start"]),
                    int(record["end"]),
                    int(record["strand"]),
                )
            )

        print(f"Prepared {len(info_rows):,} ENS_INFO rows and {len(locus_rows):,} ENS_Locus rows")
        return info_rows, locus_rows

    def ensure_go_obo(self) -> Optional[Path]:
        if self.go_obo and self.go_obo.exists():
            return self.go_obo

        if not self.download_go_obo:
            return None

        cache_dir = self.cache_dir or (Path(tempfile.gettempdir()) / "miniodp-geneinfo-cache")
        cache_dir.mkdir(parents=True, exist_ok=True)
        target = cache_dir / "go-basic.obo"
        if target.exists():
            print(f"Using cached GO ontology: {target}")
            return target

        last_error: Optional[Exception] = None
        for url in GO_BASIC_OBO_URLS:
            try:
                print(f"Downloading GO ontology from {url}")
                request = urllib.request.Request(url, headers={"User-Agent": "miniodp-geneinfo-builder/1.0"})
                with urllib.request.urlopen(request) as response, open(target, "wb") as out_handle:
                    out_handle.write(response.read())
                print(f"Saved GO ontology to {target}")
                return target
            except Exception as exc:
                last_error = exc
                logger.warning("Failed to download GO ontology from %s: %s", url, exc)

        raise RuntimeError(f"Unable to download GO ontology: {last_error}")

    def load_go_term_map(self) -> Tuple[Dict[str, Tuple[str, str, str]], Dict[str, str]]:
        go_obo_path = self.ensure_go_obo()
        if not go_obo_path:
            print("GO ontology file not provided; GO term names and namespaces will be left blank")
            return {}, {}

        term_map: Dict[str, Tuple[str, str, str]] = {}
        obsolete_map: Dict[str, str] = {}
        current: Dict[str, str] = {}
        current_alt_ids: List[str] = []

        def flush_current() -> None:
            nonlocal current_alt_ids
            term_id = current.get("id")
            if not term_id:
                return

            if current.get("is_obsolete") == "true":
                replacement = current.get("replaced_by", "").strip()
                if replacement:
                    obsolete_map[term_id] = replacement
                current_alt_ids = []
                return

            term_payload = (
                current.get("name", ""),
                extract_go_definition(current.get("def", "")),
                current.get("namespace", ""),
            )
            term_map[term_id] = term_payload
            for alt_id in current_alt_ids:
                term_map[alt_id] = term_payload
            current_alt_ids = []

        with open(go_obo_path, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if line == "[Term]":
                    flush_current()
                    current = {}
                    current_alt_ids = []
                    continue
                if not line:
                    flush_current()
                    current = {}
                    current_alt_ids = []
                    continue
                if line.startswith("["):
                    flush_current()
                    current = {}
                    current_alt_ids = []
                    continue
                if ": " in line:
                    key, value = line.split(": ", 1)
                    if key == "alt_id":
                        current_alt_ids.append(value)
                    elif key in {"id", "name", "namespace", "def", "is_obsolete", "replaced_by"}:
                        current[key] = value

        flush_current()
        print(f"Loaded {len(term_map):,} GO terms from ontology")
        print(f"Loaded {len(obsolete_map):,} obsolete GO replacements")
        return term_map, obsolete_map

    def parse_go_rows(self) -> List[Tuple[str, str, str, str, str, str]]:
        if not self.go_kegg_tsv:
            print("GO/KEGG TSV not provided; ENS_GO will be empty")
            return []

        print(f"Parsing GO annotations from {self.go_kegg_tsv}")
        go_term_map, obsolete_map = self.load_go_term_map()
        go_by_gene: Dict[str, set[str]] = defaultdict(set)

        with open_text_auto(self.go_kegg_tsv) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                gene_id = (row.get("#gene_id") or row.get("gene_id") or "").strip()
                if not gene_id:
                    continue
                raw_terms = (row.get("GO_terms") or "").strip().strip('"')
                if not raw_terms:
                    continue
                for term in raw_terms.split(","):
                    term = term.strip()
                    if term:
                        go_by_gene[gene_id].add(term)

        go_rows: List[Tuple[str, str, str, str, str, str]] = []
        for gene_id in sorted(go_by_gene):
            for go_id in sorted(go_by_gene[gene_id]):
                resolved_go_id = go_id
                seen: set[str] = set()
                while resolved_go_id in obsolete_map and resolved_go_id not in seen:
                    seen.add(resolved_go_id)
                    resolved_go_id = obsolete_map[resolved_go_id]
                go_name, go_definition, go_domain = go_term_map.get(resolved_go_id, ("", "", ""))
                go_rows.append((gene_id, go_id, go_name, go_definition, "", go_domain))

        print(f"Prepared {len(go_rows):,} ENS_GO rows across {len(go_by_gene):,} genes")
        return go_rows

    def create_schema(self, conn: sqlite3.Connection) -> None:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE ENS_INFO (
                ENS TEXT,
                GeneName TEXT,
                GeneType TEXT,
                Description TEXT
            )
            """
        )
        cursor.execute(
            """
            CREATE TABLE ENS_Locus (
                ENS TEXT,
                chromosome TEXT,
                start INTEGER,
                end INTEGER,
                strand INTEGER
            )
            """
        )
        cursor.execute(
            """
            CREATE TABLE ENS_GO (
                ENS TEXT,
                GO_Accession TEXT,
                GO_Name TEXT,
                GO_Definition TEXT,
                GO_Evidence TEXT,
                GO_Domain TEXT
            )
            """
        )

    def create_indexes(self, conn: sqlite3.Connection) -> None:
        cursor = conn.cursor()
        statements = [
            "CREATE INDEX IF NOT EXISTS idx_ens_info_ens ON ENS_INFO(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_info_genename ON ENS_INFO(GeneName)",
            "CREATE INDEX IF NOT EXISTS idx_ens_locus_ens ON ENS_Locus(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_go_ens ON ENS_GO(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_go_accession ON ENS_GO(GO_Accession)",
            "CREATE INDEX IF NOT EXISTS idx_ens_go_name ON ENS_GO(GO_Name)",
        ]
        for statement in statements:
            cursor.execute(statement)

    def build_database(self) -> Tuple[int, int]:
        info_rows, locus_rows = self.parse_gtf()
        go_rows = self.parse_go_rows()

        if self.geneinfo_db.exists():
            self.geneinfo_db.unlink()

        with sqlite3.connect(self.geneinfo_db) as conn:
            self.create_schema(conn)
            cursor = conn.cursor()
            cursor.executemany("INSERT INTO ENS_INFO VALUES (?, ?, ?, ?)", info_rows)
            cursor.executemany("INSERT INTO ENS_Locus VALUES (?, ?, ?, ?, ?)", locus_rows)
            if go_rows:
                cursor.executemany("INSERT INTO ENS_GO VALUES (?, ?, ?, ?, ?, ?)", go_rows)
            conn.commit()
            self.create_indexes(conn)
            cursor.execute("ANALYZE")
            conn.commit()

            gene_count = cursor.execute("SELECT COUNT(*) FROM ENS_INFO").fetchone()[0]
            go_count = cursor.execute("SELECT COUNT(*) FROM ENS_GO").fetchone()[0]
            if gene_count <= 0:
                raise RuntimeError("Database built but ENS_INFO is empty")

        print(f"Built geneinfo.db with {gene_count:,} genes")
        return gene_count, go_count

    def generate_toml_config(self) -> None:
        print("Generating geneinfo.toml configuration")

        has_descriptions = False
        has_go_rows = False
        has_go_names = False
        with sqlite3.connect(self.geneinfo_db) as conn:
            cursor = conn.cursor()
            has_descriptions = bool(
                cursor.execute(
                    "SELECT 1 FROM ENS_INFO WHERE Description IS NOT NULL AND TRIM(Description) != '' LIMIT 1"
                ).fetchone()
            )
            has_go_rows = bool(cursor.execute("SELECT 1 FROM ENS_GO LIMIT 1").fetchone())
            has_go_names = bool(
                cursor.execute(
                    "SELECT 1 FROM ENS_GO WHERE GO_Name IS NOT NULL AND TRIM(GO_Name) != '' LIMIT 1"
                ).fetchone()
            )

        lines = [
            "# Geneinfo database search configuration",
            f"# Auto-generated for species: {self.species_label}",
            f"# Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "[[search_fields]]",
            'label = "Gene Name"',
            'value = "gene_name"',
            'db_table = "ENS_INFO"',
            'db_column = "GeneName"',
            'match_type = "like"',
            "",
            "[[search_fields]]",
            'label = "Gene ID"',
            'value = "ens_id"',
            'db_table = "ENS_INFO"',
            'db_column = "ENS"',
            'match_type = "exact"',
            "",
        ]

        if has_descriptions:
            lines.extend(
                [
                    "[[search_fields]]",
                    'label = "Description"',
                    'value = "description"',
                    'db_table = "ENS_INFO"',
                    'db_column = "Description"',
                    'match_type = "like"',
                    "",
                ]
            )

        if has_go_rows:
            lines.extend(
                [
                    "[[search_fields]]",
                    'label = "GO Accession"',
                    'value = "go_id"',
                    'db_table = "ENS_GO"',
                    'db_column = "GO_Accession"',
                    'match_type = "exact"',
                    "",
                ]
            )

        if has_go_names:
            lines.extend(
                [
                    "[[search_fields]]",
                    'label = "GO Term"',
                    'value = "go_term"',
                    'db_table = "ENS_GO"',
                    'db_column = "GO_Name"',
                    'match_type = "like"',
                    "",
                ]
            )

        self.geneinfo_toml.write_text("\n".join(lines), encoding="utf-8")

    def run(self) -> None:
        self.validate_inputs()
        self.build_database()
        self.generate_toml_config()
        print("Output files:")
        print(f"  - {self.geneinfo_db}")
        print(f"  - {self.geneinfo_toml}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build geneinfo.db from local annotation files",
        epilog=(
            "Example: python dash/scripts/build_geneinfo_from_local.py "
            "--gtf annotation.gtf --annotation-tsv annotation.tsv --go-kegg-tsv gene2go.tsv "
            "--download-go-obo -o ./build/geneinfo"
        ),
    )
    parser.add_argument("--gtf", required=True, help="Input annotation GTF/GFF path")
    parser.add_argument(
        "--annotation-tsv",
        help="Optional gene-level annotation TSV/TSV.GZ/ZIP with gene symbol and description columns",
    )
    parser.add_argument("--go-kegg-tsv", help="Optional gene-to-GO/KEGG TSV path")
    parser.add_argument("--go-obo", help="Optional local go-basic.obo path")
    parser.add_argument(
        "--cache-dir",
        help="Optional cache directory for downloaded auxiliary files such as go-basic.obo",
    )
    parser.add_argument(
        "--download-go-obo",
        action="store_true",
        help="Download go-basic.obo automatically if --go-obo is not provided",
    )
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory for geneinfo.db and geneinfo.toml")
    parser.add_argument("--species-label", default="local_species", help="Species label written into geneinfo.toml")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    builder = LocalGeneinfoBuilder(
        gtf_file=Path(args.gtf),
        annotation_tsv=Path(args.annotation_tsv) if args.annotation_tsv else None,
        go_kegg_tsv=Path(args.go_kegg_tsv) if args.go_kegg_tsv else None,
        go_obo=Path(args.go_obo) if args.go_obo else None,
        cache_dir=Path(args.cache_dir) if args.cache_dir else None,
        download_go_obo=args.download_go_obo,
        output_dir=Path(args.output_dir),
        species_label=args.species_label,
    )

    try:
        builder.run()
    except Exception as exc:
        logger.error("Build failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
