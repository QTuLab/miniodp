#!/usr/bin/env python3
"""Summarize bulk run statistics from per-species CSV files.

Each CSV represents one species and contains run-level metadata. The script
counts non-empty Run rows, sums Spots, optionally reconstructs Bulk bases from
formal evidence tables, and updates species_display.toml while preserving
comments and order.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from dataclasses import dataclass, field
from datetime import datetime
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable

try:
    import tomlkit  # type: ignore
except ImportError as exc:  # pragma: no cover - dependency is expected locally
    raise RuntimeError(
        "'tomlkit' is required to preserve species_display.toml formatting."
    ) from exc

from logging_utils import configure_script_logger, create_print_proxy  # type: ignore


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

_SCRIPT_LOGGER = configure_script_logger(__name__)
print = create_print_proxy(_SCRIPT_LOGGER)


REQUIRED_COLUMNS = ("Study", "Sample", "Library Strategy", "Run", "Spots")
PUBLIC_BASES_REQUIRED_COLUMNS = ("species", "public_summary_tsv")
INTERNAL_BASES_REQUIRED_COLUMNS = (
    "species",
    "run",
    "spots",
    "bases",
    "evidence_source",
    "evidence_detail",
    "notes",
)
RUN_ID_RE = re.compile(r"^[SED]RR\d+$")
PUBLIC_RUN_PREFIXES = ("SRR", "ERR", "DRR")


@dataclass
class BulkStat:
    """Aggregated bulk run statistics for one species."""

    species_key: str
    source_file: Path
    run_count: int
    spots_sum: int
    public_bases: int | None = None
    internal_bases: int = 0
    missing_internal_base_rows: list[str] = field(default_factory=list)
    current_bulk_datasets: int | None = None
    current_reads: int | None = None
    current_bulk_bases: int | None = None
    run_ids: set[str] = field(default_factory=set)
    status: str = "pending"
    notes: list[str] = field(default_factory=list)

    @property
    def new_bulk_datasets(self) -> int:
        return self.run_count

    @property
    def new_reads(self) -> int:
        return self.spots_sum

    @property
    def new_bulk_bases(self) -> int | None:
        if self.public_bases is None:
            return None
        if self.missing_internal_base_rows:
            return None
        return self.public_bases + self.internal_bases


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize bulk runs and spots from per-species CSV files.",
    )
    parser.add_argument(
        "--source-dir",
        required=True,
        help="Directory containing one {species_key}.csv file per species.",
    )
    parser.add_argument(
        "--species-display",
        default="hugo/data/species_display.toml",
        help="Hugo species_display.toml path to inspect or update.",
    )
    parser.add_argument(
        "--species",
        help="Optional species key to update, for example medaka.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned updates without writing species_display.toml.",
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Create a .bak copy before writing species_display.toml.",
    )
    parser.add_argument(
        "--output-csv",
        help="Optional audit CSV path.",
    )
    parser.add_argument(
        "--public-bases-table",
        default=None,
        help="Optional TSV mapping species to audited public bases summary files.",
    )
    parser.add_argument(
        "--internal-bases-table",
        default=None,
        help="Optional TSV with run-level bases evidence for internal nonpublic bulk runs.",
    )
    parser.add_argument(
        "--strict-run-id",
        action="store_true",
        help="Report non-standard run IDs such as private GEO sample IDs.",
    )
    return parser.parse_args()


def normalize_header(value: object) -> str:
    return str(value).strip() if value is not None else ""


def normalize_cell(value: object) -> str:
    return str(value).strip() if value is not None else ""


def normalize_int(value: object) -> int | None:
    text = normalize_cell(value)
    if not text:
        return None
    cleaned = text.replace(",", "").replace(" ", "")
    try:
        number = Decimal(cleaned)
    except InvalidOperation as exc:
        raise ValueError(f"Value is not numeric: {value!r}") from exc
    if number != number.to_integral_value():
        raise ValueError(f"Value is not an integer: {value!r}")
    return int(number)


def parse_spots(value: object, source_file: Path, row_number: int) -> int:
    """Parse Spots as an integer while accepting comma-formatted strings."""

    text = normalize_cell(value)
    if not text:
        raise ValueError(f"{source_file} row {row_number}: Spots is missing")

    cleaned = text.replace(",", "").replace(" ", "")
    try:
        number = Decimal(cleaned)
    except InvalidOperation as exc:
        raise ValueError(
            f"{source_file} row {row_number}: Spots is not numeric: {value!r}"
        ) from exc

    if number != number.to_integral_value():
        raise ValueError(
            f"{source_file} row {row_number}: Spots is not an integer: {value!r}"
        )
    return int(number)


def format_row_list(rows: Iterable[int]) -> str:
    return ", ".join(str(row) for row in rows)


def find_csv_files(source_dir: Path, selected_species: str | None) -> list[Path]:
    if not source_dir.exists():
        raise FileNotFoundError(f"Source directory not found: {source_dir}")
    if not source_dir.is_dir():
        raise NotADirectoryError(f"Source path is not a directory: {source_dir}")

    if selected_species:
        path = source_dir / f"{selected_species}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Species CSV not found: {path}")
        return [path]

    csv_files = sorted(source_dir.glob("*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No species CSV files found in {source_dir}")
    return csv_files


def validate_columns(source_file: Path, fieldnames: list[str] | None) -> None:
    headers = [normalize_header(value) for value in fieldnames or []]
    missing = [column for column in REQUIRED_COLUMNS if column not in headers]
    if missing:
        raise ValueError(
            f"{source_file}: required columns missing: {', '.join(missing)}"
        )


def validate_named_columns(
    source_file: Path,
    fieldnames: list[str] | None,
    required_columns: tuple[str, ...],
) -> None:
    headers = [normalize_header(value) for value in fieldnames or []]
    missing = [column for column in required_columns if column not in headers]
    if missing:
        raise ValueError(
            f"{source_file}: required columns missing: {', '.join(missing)}"
        )


def row_has_content(row: dict[str, str | None]) -> bool:
    return any(normalize_cell(value) for value in row.values())


def summarize_csv(source_file: Path, strict_run_id: bool) -> BulkStat:
    species_key = source_file.stem
    run_seen: dict[str, int] = {}
    invalid_runs: list[str] = []
    empty_run_rows: list[int] = []

    run_count = 0
    spots_sum = 0
    run_ids: set[str] = set()

    with open(source_file, "r", newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile)
        validate_columns(source_file, reader.fieldnames)

        for row_number, row in enumerate(reader, start=2):
            if not row_has_content(row):
                continue

            run = normalize_cell(row.get("Run"))
            if not run:
                empty_run_rows.append(row_number)
                continue

            if run in run_seen:
                raise ValueError(
                    f"{source_file} row {row_number}: duplicate Run {run!r}; "
                    f"first seen at row {run_seen[run]}"
                )
            run_seen[run] = row_number

            if strict_run_id and not RUN_ID_RE.match(run):
                invalid_runs.append(f"{run} (row {row_number})")

            run_ids.add(run)

            spots = parse_spots(row.get("Spots"), source_file, row_number)
            run_count += 1
            spots_sum += spots

    if empty_run_rows:
        raise ValueError(
            f"{source_file}: non-empty rows have empty Run values: "
            f"{format_row_list(empty_run_rows)}"
        )

    stat = BulkStat(
        species_key=species_key,
        source_file=source_file,
        run_count=run_count,
        spots_sum=spots_sum,
        run_ids=run_ids,
    )

    if strict_run_id and invalid_runs:
        stat.notes.append("non_standard_run_id=" + "; ".join(invalid_runs))

    return stat


def load_stats(
    source_dir: Path,
    selected_species: str | None,
    strict_run_id: bool,
) -> list[BulkStat]:
    return [
        summarize_csv(source_file, strict_run_id)
        for source_file in find_csv_files(source_dir, selected_species)
    ]


def load_public_bases_map(path: Path | None) -> dict[str, tuple[int | None, str, list[str]]]:
    if path is None or not path.exists():
        return {}

    public_bases_map: dict[str, tuple[int | None, str, list[str]]] = {}
    with open(path, "r", newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        validate_named_columns(path, reader.fieldnames, PUBLIC_BASES_REQUIRED_COLUMNS)
        for row_number, row in enumerate(reader, start=2):
            species = normalize_cell(row.get("species"))
            if not species:
                continue
            summary_ref = normalize_cell(row.get("public_summary_tsv"))
            if not summary_ref:
                raise ValueError(f"{path} row {row_number}: public_summary_tsv is missing")
            summary_path = Path(summary_ref)
            if not summary_path.is_absolute():
                summary_path = (Path.cwd() / summary_path).resolve()
            notes: list[str] = []
            if not summary_path.exists():
                notes.append(f"public_bases_summary_missing={summary_path}")
                public_bases_map[species] = (None, str(summary_path), notes)
                continue
            with open(summary_path, "r", newline="", encoding="utf-8-sig") as handle:
                summary_reader = csv.DictReader(handle, delimiter="\t")
                matching_row = None
                for summary_row in summary_reader:
                    if normalize_cell(summary_row.get("species")) == species:
                        matching_row = summary_row
                        break
            if matching_row is None:
                notes.append(f"public_bases_summary_missing_species={summary_path}")
                public_bases_map[species] = (None, str(summary_path), notes)
                continue
            audit_status = normalize_cell(matching_row.get("audit_status"))
            if audit_status and audit_status != "pass":
                notes.append(f"public_bases_audit_status={audit_status}")
                public_bases_map[species] = (None, str(summary_path), notes)
                continue
            bases_value = matching_row.get("expected_bases")
            try:
                public_bases = normalize_int(bases_value)
            except ValueError as exc:
                raise ValueError(
                    f"{summary_path}: expected_bases is invalid: {bases_value!r}"
                ) from exc
            if public_bases is None:
                notes.append("public_bases_missing")
            public_bases_map[species] = (public_bases, str(summary_path), notes)
    return public_bases_map


def load_internal_bases_map(path: Path | None) -> dict[str, dict[str, object]]:
    if path is None or not path.exists():
        return {}

    grouped: dict[str, dict[str, object]] = {}
    with open(path, "r", newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        validate_named_columns(path, reader.fieldnames, INTERNAL_BASES_REQUIRED_COLUMNS)
        for row_number, row in enumerate(reader, start=2):
            species = normalize_cell(row.get("species"))
            if not species:
                continue
            run = normalize_cell(row.get("run"))
            if not run:
                raise ValueError(f"{path} row {row_number}: run is missing")
            grouped.setdefault(
                species,
                {"bases_sum": 0, "missing_rows": [], "notes": []},
            )
            species_entry = grouped[species]
            raw_bases = row.get("bases")
            if not normalize_cell(raw_bases):
                species_entry["missing_rows"].append(f"{run} (row {row_number})")
                continue
            try:
                bases = normalize_int(raw_bases)
            except ValueError as exc:
                raise ValueError(
                    f"{path} row {row_number}: bases is invalid: {raw_bases!r}"
                ) from exc
            if bases is None:
                species_entry["missing_rows"].append(f"{run} (row {row_number})")
                continue
            species_entry["bases_sum"] += bases
            species_entry.setdefault("run_bases", {})[run] = bases
    return grouped


def load_toml(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"species_display.toml not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        return tomlkit.parse(fh.read())


def attach_current_values(stats: Iterable[BulkStat], toml_data) -> None:
    for stat in stats:
        if stat.species_key not in toml_data:
            raise ValueError(
                f"{stat.source_file}: species key {stat.species_key!r} "
                "not found in species_display.toml"
            )

        section = toml_data[stat.species_key]
        stat.current_bulk_datasets = int(section.get("bulk_datasets", 0))
        stat.current_reads = int(section.get("reads", 0))
        bulk_bases_value = section.get("bulk_bases")
        if bulk_bases_value is not None:
            stat.current_bulk_bases = int(bulk_bases_value)
        if (
            stat.current_bulk_datasets == stat.new_bulk_datasets
            and stat.current_reads == stat.new_reads
            and (
                stat.new_bulk_bases is None
                or stat.current_bulk_bases == stat.new_bulk_bases
            )
        ):
            stat.status = "unchanged"
        else:
            stat.status = "changed"


def attach_bulk_bases_evidence(
    stats: Iterable[BulkStat],
    public_bases_map: dict[str, tuple[int | None, str, list[str]]],
    internal_bases_map: dict[str, dict[str, object]],
) -> None:
    for stat in stats:
        public_bases_entry = public_bases_map.get(stat.species_key)
        if public_bases_entry is not None:
            public_bases, summary_path, notes = public_bases_entry
            stat.public_bases = public_bases
            if summary_path:
                stat.notes.append(f"public_bases_summary={summary_path}")
            stat.notes.extend(notes)

        internal_entry = internal_bases_map.get(stat.species_key)
        if internal_entry is not None:
            run_bases = {
                str(run): int(value)
                for run, value in dict(internal_entry.get("run_bases", {})).items()
            }
            stat.internal_bases = sum(value for run, value in run_bases.items() if run in stat.run_ids)
            missing_rows = []
            for item in list(internal_entry["missing_rows"]):
                run_id = str(item).split(" ", 1)[0]
                if run_id in stat.run_ids:
                    missing_rows.append(item)
            stat.missing_internal_base_rows = missing_rows
            if stat.missing_internal_base_rows:
                stat.notes.append(
                    "missing_internal_bases="
                    + "; ".join(stat.missing_internal_base_rows)
                )


def print_summary(stats: Iterable[BulkStat]) -> None:
    rows = list(stats)
    print("species          bulk_datasets old -> new        reads old -> new        bulk_bases old -> new")
    for stat in rows:
        old_bulk = "NA" if stat.current_bulk_datasets is None else str(stat.current_bulk_datasets)
        old_reads = "NA" if stat.current_reads is None else str(stat.current_reads)
        old_bases = "NA" if stat.current_bulk_bases is None else str(stat.current_bulk_bases)
        new_bases = "pending" if stat.new_bulk_bases is None else str(stat.new_bulk_bases)
        print(
            f"{stat.species_key:<16} "
            f"{old_bulk:>8} -> {stat.new_bulk_datasets:<8} "
            f"{old_reads:>18} -> {stat.new_reads:<18} "
            f"{old_bases:>18} -> {new_bases}"
        )
        if stat.notes:
            print(f"  notes: {' | '.join(stat.notes)}")


def write_audit_csv(path: Path, stats: Iterable[BulkStat]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cwd = Path.cwd()
    with open(path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "species_key",
                "source_file",
                "run_count",
                "spots_sum",
                "public_bases",
                "internal_bases",
                "missing_internal_base_rows",
                "current_bulk_datasets",
                "current_reads",
                "current_bulk_bases",
                "new_bulk_datasets",
                "new_reads",
                "new_bulk_bases",
                "status",
                "notes",
            ]
        )
        for stat in stats:
            try:
                source_file = str(stat.source_file.relative_to(cwd))
            except ValueError:
                source_file = str(stat.source_file)
            writer.writerow(
                [
                    stat.species_key,
                    source_file,
                    stat.run_count,
                    stat.spots_sum,
                    stat.public_bases if stat.public_bases is not None else "",
                    stat.internal_bases,
                    " | ".join(stat.missing_internal_base_rows),
                    stat.current_bulk_datasets if stat.current_bulk_datasets is not None else "",
                    stat.current_reads if stat.current_reads is not None else "",
                    stat.current_bulk_bases if stat.current_bulk_bases is not None else "",
                    stat.new_bulk_datasets,
                    stat.new_reads,
                    stat.new_bulk_bases if stat.new_bulk_bases is not None else "",
                    stat.status,
                    " | ".join(stat.notes),
                ]
            )
    print(f"Audit CSV written to {path}")


def update_species_display(
    path: Path,
    toml_data,
    stats: Iterable[BulkStat],
    dry_run: bool,
    backup: bool,
) -> None:
    writable_stats = [stat for stat in stats if stat.status != "missing_species"]
    if dry_run:
        print("Dry run: no changes written to species_display.toml")
        return

    if not writable_stats:
        raise RuntimeError("No matching species sections found to update")

    if backup:
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        backup_dir = path.parent.parent / f"{path.parent.name}_backups"
        backup_dir.mkdir(parents=True, exist_ok=True)
        backup_path = backup_dir / f"{path.name}.{timestamp}.bak"
        shutil.copy2(path, backup_path)
        print(f"Backup written to {backup_path}")

    for stat in writable_stats:
        section = toml_data[stat.species_key]
        section["bulk_datasets"] = tomlkit.integer(stat.new_bulk_datasets)
        section["reads"] = tomlkit.integer(stat.new_reads)
        if stat.new_bulk_bases is not None:
            section["bulk_bases"] = tomlkit.integer(stat.new_bulk_bases)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(tomlkit.dumps(toml_data))
    print(f"species_display.toml updated: {path}")


def main() -> None:
    args = parse_args()
    source_dir = Path(args.source_dir).resolve()
    species_display = Path(args.species_display).resolve()
    output_csv = Path(args.output_csv).resolve() if args.output_csv else None
    public_bases_table = Path(args.public_bases_table).resolve() if args.public_bases_table else None
    internal_bases_table = Path(args.internal_bases_table).resolve() if args.internal_bases_table else None

    stats = load_stats(source_dir, args.species, args.strict_run_id)
    public_bases_map = load_public_bases_map(public_bases_table)
    internal_bases_map = load_internal_bases_map(internal_bases_table)
    attach_bulk_bases_evidence(stats, public_bases_map, internal_bases_map)
    toml_data = load_toml(species_display)
    attach_current_values(stats, toml_data)
    print_summary(stats)

    if output_csv:
        write_audit_csv(output_csv, stats)

    update_species_display(
        path=species_display,
        toml_data=toml_data,
        stats=stats,
        dry_run=args.dry_run,
        backup=args.backup,
    )


if __name__ == "__main__":
    main()
