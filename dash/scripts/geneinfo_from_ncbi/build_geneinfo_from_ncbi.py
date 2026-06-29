#!/usr/bin/env python3
"""
Build geneinfo.db from local NCBI RefSeq files.

This script creates a miniodp-compatible geneinfo database using:
- NCBI feature table (.txt.gz)
- NCBI RNA FASTA (.fna)
- Optional NCBI GO GAF (.gaf.gz)
- Optional GO OBO ontology file for GO term names/definitions
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sqlite3
import sys
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

GO_BASIC_OBO_URLS = [
    "https://current.geneontology.org/ontology/go-basic.obo",
    "https://purl.obolibrary.org/obo/go/go-basic.obo",
]


def open_text_auto(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def extract_go_definition(raw_value: str) -> str:
    raw_value = raw_value.strip()
    if raw_value.startswith('"') and '"' in raw_value[1:]:
        return raw_value.split('"', 2)[1]
    return raw_value


class NCBIGeneinfoBuilder:
    def __init__(
        self,
        *,
        feature_table: Path,
        rna_fasta: Path,
        output_dir: Path,
        go_gaf: Optional[Path] = None,
        go_obo: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        download_go_obo: bool = False,
        species_label: str = "killifish",
    ) -> None:
        self.feature_table = feature_table
        self.rna_fasta = rna_fasta
        self.go_gaf = go_gaf
        self.go_obo = go_obo
        self.cache_dir = cache_dir
        self.download_go_obo = download_go_obo
        self.species_label = species_label
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.geneinfo_db = self.output_dir / "geneinfo.db"
        self.geneinfo_toml = self.output_dir / "geneinfo.toml"

    def validate_inputs(self) -> None:
        required = {
            "feature_table": self.feature_table,
            "rna_fasta": self.rna_fasta,
        }
        for label, path in required.items():
            if not path.exists():
                raise FileNotFoundError(f"Missing required input {label}: {path}")
        if self.go_gaf and not self.go_gaf.exists():
            raise FileNotFoundError(f"GO GAF file not found: {self.go_gaf}")
        if self.go_obo and not self.go_obo.exists():
            raise FileNotFoundError(f"GO OBO file not found: {self.go_obo}")

    def iter_feature_rows(self) -> Iterator[Dict[str, str]]:
        with open_text_auto(self.feature_table) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if not row:
                    continue
                normalized = dict(row)
                if "# feature" in normalized:
                    normalized["feature"] = normalized.pop("# feature")
                yield normalized

    def build_gene_tables(self) -> Tuple[List[Tuple[str, str, str, str]], List[Tuple[str, str, int, int, int]]]:
        info_rows: List[Tuple[str, str, str, str]] = []
        locus_rows: List[Tuple[str, str, int, int, int]] = []
        seen_gene_ids = set()

        for row in self.iter_feature_rows():
            if row.get("feature") != "gene":
                continue

            gene_id = (row.get("GeneID") or "").strip()
            if not gene_id or gene_id in seen_gene_ids:
                continue
            seen_gene_ids.add(gene_id)

            gene_name = (row.get("symbol") or "").strip() or gene_id
            gene_type = (row.get("class") or "").strip() or "unknown"
            description = (row.get("name") or "").strip() or gene_name
            chromosome = (row.get("chromosome") or "").strip() or (row.get("genomic_accession") or "").strip()
            start = int(row.get("start") or 0)
            end = int(row.get("end") or 0)
            strand_symbol = (row.get("strand") or "").strip()
            strand = 1 if strand_symbol == "+" else -1 if strand_symbol == "-" else 0

            info_rows.append((gene_id, gene_name, gene_type, description))
            locus_rows.append((gene_id, chromosome, start, end, strand))

        logger.info("Prepared %s ENS_INFO rows and %s ENS_Locus rows", len(info_rows), len(locus_rows))
        return info_rows, locus_rows

    def build_transcript_index(self) -> Dict[str, str]:
        transcript_to_gene: Dict[str, str] = {}

        for row in self.iter_feature_rows():
            feature = row.get("feature", "")
            if feature not in {"mRNA", "ncRNA", "rRNA", "tRNA", "transcript"}:
                continue

            gene_id = (row.get("GeneID") or "").strip()
            transcript_id = (row.get("product_accession") or "").strip()
            if gene_id and transcript_id and transcript_id not in transcript_to_gene:
                transcript_to_gene[transcript_id] = gene_id

        logger.info("Prepared transcript index for %s transcript accessions", len(transcript_to_gene))
        return transcript_to_gene

    def iter_fasta(self, path: Path) -> Iterator[Tuple[str, str]]:
        header: Optional[str] = None
        seq_parts: List[str] = []
        with open_text_auto(path) as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_parts)
                    header = line[1:]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if header is not None:
                yield header, "".join(seq_parts)

    def build_seq_rows(self, transcript_to_gene: Dict[str, str]) -> List[Tuple[str, str, str]]:
        seq_rows: List[Tuple[str, str, str]] = []
        missing = 0

        for header, sequence in self.iter_fasta(self.rna_fasta):
            transcript_id = header.split(None, 1)[0]
            gene_id = transcript_to_gene.get(transcript_id)
            if not gene_id:
                missing += 1
                continue
            seq_rows.append((gene_id, transcript_id, sequence))

        logger.info(
            "Prepared %s ENS_Seq rows from RNA FASTA (%s headers had no GeneID mapping)",
            len(seq_rows),
            missing,
        )
        return seq_rows

    def ensure_go_obo(self) -> Optional[Path]:
        if self.go_obo and self.go_obo.exists():
            return self.go_obo

        if not self.download_go_obo:
            return None

        cache_dir = self.cache_dir or (Path(tempfile.gettempdir()) / "miniodp-geneinfo-cache")
        cache_dir.mkdir(parents=True, exist_ok=True)
        target = cache_dir / "go-basic.obo"
        if target.exists():
            logger.info("Using cached GO ontology: %s", target)
            return target

        last_error: Optional[Exception] = None
        for url in GO_BASIC_OBO_URLS:
            try:
                logger.info("Downloading GO ontology from %s", url)
                request = urllib.request.Request(url, headers={"User-Agent": "miniodp-geneinfo-builder/1.0"})
                with urllib.request.urlopen(request) as response, open(target, "wb") as out_handle:
                    out_handle.write(response.read())
                logger.info("Saved GO ontology to %s", target)
                return target
            except Exception as exc:
                last_error = exc
                logger.warning("Failed to download GO ontology from %s: %s", url, exc)

        raise RuntimeError(f"Unable to download GO ontology: {last_error}")

    def load_go_term_map(self) -> Dict[str, Tuple[str, str, str]]:
        go_obo_path = self.ensure_go_obo()
        if not go_obo_path:
            logger.warning("No GO ontology file available; GO term names/definitions will be left blank")
            return {}

        term_map: Dict[str, Tuple[str, str, str]] = {}
        current: Dict[str, str] = {}

        def flush_current() -> None:
            term_id = current.get("id")
            if not term_id or current.get("is_obsolete") == "true":
                return
            name = current.get("name", "")
            definition = extract_go_definition(current.get("def", ""))
            namespace = current.get("namespace", "")
            term_map[term_id] = (name, definition, namespace)

        with open(go_obo_path, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if line == "[Term]":
                    flush_current()
                    current = {}
                    continue
                if not line:
                    flush_current()
                    current = {}
                    continue
                if line.startswith("["):
                    flush_current()
                    current = {}
                    continue
                if ": " in line:
                    key, value = line.split(": ", 1)
                    if key in {"id", "name", "namespace", "def", "is_obsolete"}:
                        current[key] = value

        logger.info("Loaded %s GO terms from ontology file", len(term_map))
        return term_map

    def build_go_rows(self) -> List[Tuple[str, str, str, str, str, str]]:
        if not self.go_gaf:
            logger.info("No GO GAF provided; ENS_GO table will be empty")
            return []

        go_term_map = self.load_go_term_map()
        domain_map = {
            "P": "biological_process",
            "F": "molecular_function",
            "C": "cellular_component",
        }

        rows: List[Tuple[str, str, str, str, str, str]] = []
        seen = set()

        with open_text_auto(self.go_gaf) as handle:
            reader = csv.reader(handle, delimiter="\t")
            for fields in reader:
                if not fields or fields[0].startswith("!"):
                    continue
                if len(fields) < 9:
                    continue

                gene_id = fields[1].strip()
                go_id = fields[4].strip()
                evidence = fields[6].strip() or "IEA"
                aspect = fields[8].strip()
                if not gene_id or not go_id:
                    continue

                term_name, definition, namespace = go_term_map.get(go_id, ("", "", ""))
                go_domain = namespace or domain_map.get(aspect, "")
                row = (gene_id, go_id, term_name, definition, evidence, go_domain)
                if row in seen:
                    continue
                seen.add(row)
                rows.append(row)

        logger.info("Prepared %s ENS_GO rows", len(rows))
        return rows

    def create_schema(self, conn: sqlite3.Connection) -> None:
        cursor = conn.cursor()
        cursor.executescript(
            """
            CREATE TABLE ENS_INFO (
                ENS TEXT,
                GeneName TEXT,
                GeneType TEXT,
                Description TEXT
            );

            CREATE TABLE ENS_Locus (
                ENS TEXT,
                chromosome TEXT,
                start INTEGER,
                end INTEGER,
                strand INTEGER
            );

            CREATE TABLE ENS_GO (
                ENS TEXT,
                GO_Accession TEXT,
                GO_Name TEXT,
                GO_Definition TEXT,
                GO_Evidence TEXT,
                GO_Domain TEXT
            );

            CREATE TABLE ENS_Seq (
                ENS TEXT,
                ENStrpt TEXT,
                nseq TEXT
            );
            """
        )
        conn.commit()

    def create_indexes(self, conn: sqlite3.Connection) -> None:
        cursor = conn.cursor()
        for sql in [
            "CREATE INDEX IF NOT EXISTS idx_ens_info_ens ON ENS_INFO(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_info_genename ON ENS_INFO(GeneName)",
            "CREATE INDEX IF NOT EXISTS idx_ens_locus_ens ON ENS_Locus(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_go_ens ON ENS_GO(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_go_accession ON ENS_GO(GO_Accession)",
            "CREATE INDEX IF NOT EXISTS idx_ens_go_name ON ENS_GO(GO_Name)",
            "CREATE INDEX IF NOT EXISTS idx_ens_seq_ens ON ENS_Seq(ENS)",
        ]:
            cursor.execute(sql)
        conn.commit()

    def build_database(self) -> bool:
        info_rows, locus_rows = self.build_gene_tables()
        transcript_index = self.build_transcript_index()
        seq_rows = self.build_seq_rows(transcript_index)
        go_rows = self.build_go_rows()

        if self.geneinfo_db.exists():
            self.geneinfo_db.unlink()

        with sqlite3.connect(self.geneinfo_db) as conn:
            self.create_schema(conn)
            cursor = conn.cursor()
            cursor.executemany("INSERT INTO ENS_INFO VALUES (?, ?, ?, ?)", info_rows)
            cursor.executemany("INSERT INTO ENS_Locus VALUES (?, ?, ?, ?, ?)", locus_rows)
            if go_rows:
                cursor.executemany("INSERT INTO ENS_GO VALUES (?, ?, ?, ?, ?, ?)", go_rows)
            if seq_rows:
                cursor.executemany("INSERT INTO ENS_Seq VALUES (?, ?, ?)", seq_rows)
            conn.commit()
            self.create_indexes(conn)
            cursor.execute("ANALYZE")
            conn.commit()

            gene_count = cursor.execute("SELECT COUNT(*) FROM ENS_INFO").fetchone()[0]
            if gene_count <= 0:
                raise RuntimeError("Database built but ENS_INFO is empty")

            logger.info("Built geneinfo.db with %s genes", gene_count)
        return True

    def generate_toml_config(self) -> bool:
        logger.info("Generating geneinfo.toml configuration")

        go_term_available = bool(self.go_gaf and self.ensure_go_obo())
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
            "[[search_fields]]",
            'label = "Description"',
            'value = "description"',
            'db_table = "ENS_INFO"',
            'db_column = "Description"',
            'match_type = "like"',
            "",
        ]

        if self.go_gaf:
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
        if go_term_available:
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
        return True

    def run(self) -> bool:
        self.validate_inputs()
        self.build_database()
        self.generate_toml_config()
        logger.info("Output files:")
        logger.info("  - %s", self.geneinfo_db)
        logger.info("  - %s", self.geneinfo_toml)
        return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build geneinfo.db from local NCBI RefSeq files",
        epilog=(
            "Example: python build_geneinfo_from_ncbi.py "
            "--feature-table feature_table.txt.gz --rna-fasta rna.fna "
            "--go-gaf gene_ontology.gaf.gz -o ./publish/dash"
        ),
    )
    parser.add_argument("--feature-table", required=True, help="NCBI feature_table.txt(.gz) path")
    parser.add_argument("--rna-fasta", required=True, help="NCBI rna.fna(.gz) path")
    parser.add_argument("--go-gaf", help="Optional NCBI GO GAF (.gaf.gz) path")
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
    parser.add_argument("--species-label", default="killifish", help="Species label written into geneinfo.toml")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    builder = NCBIGeneinfoBuilder(
        feature_table=Path(args.feature_table),
        rna_fasta=Path(args.rna_fasta),
        go_gaf=Path(args.go_gaf) if args.go_gaf else None,
        go_obo=Path(args.go_obo) if args.go_obo else None,
        cache_dir=Path(args.cache_dir) if args.cache_dir else None,
        download_go_obo=args.download_go_obo,
        output_dir=Path(args.output_dir),
        species_label=args.species_label,
    )

    try:
        success = builder.run()
    except Exception as exc:
        logger.error("Build failed: %s", exc)
        return 1
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
