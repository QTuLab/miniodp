"""Summarize single-cell dataset statistics for a species.

This script scans the `scRNA` and `scATAC` directories within a species data
folder, counts valid samples (those containing a `{sample}_metadata.json`),
aggregates the `n_cells` field, writes a per-sample CSV, prints assay-level
summaries, and updates `dash/data/species_display.toml` with the latest
single-cell totals.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:  # pragma: no cover - optional dependency for comment preservation
    import tomlkit  # type: ignore
except ImportError:  # pragma: no cover - fall back to toml
    tomlkit = None

if tomlkit is None:
    try:
        import toml  # type: ignore
    except ImportError as exc:  # pragma: no cover - dependency declared in requirements
        raise RuntimeError("Either 'tomlkit' or 'toml' is required to run this script.") from exc
else:
    toml = None  # type: ignore


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from logging_utils import configure_script_logger, create_print_proxy  # type: ignore

_SCRIPT_LOGGER = configure_script_logger(__name__)
print = create_print_proxy(_SCRIPT_LOGGER)


@dataclass
class SampleStat:
    """Record summarizing a single sample."""

    species: str
    assay: str
    study: str
    sample: str
    cells: int


def load_toml_file(path: Path):
    """Load a TOML document, preserving comments if tomlkit is available."""

    if tomlkit is not None:
        with open(path, "r", encoding="utf-8") as fh:
            return tomlkit.parse(fh.read())

    with open(path, "r", encoding="utf-8") as fh:
        return toml.load(fh)


def dump_toml_file(data, path: Path) -> None:
    """Write TOML data, preserving comments when using tomlkit."""

    if tomlkit is not None:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(tomlkit.dumps(data))
        return

    with open(path, "w", encoding="utf-8") as fh:
        toml.dump(data, fh)


def make_toml_int(value: int):
    """Return an integer value suitable for tomlkit or plain toml."""

    if tomlkit is not None:
        return tomlkit.integer(value)
    return value


def replace_or_append(section: str, key: str, value: int) -> str:
    """Replace the first occurrence of key=value in the section, or append if missing."""

    pattern = re.compile(rf"^({re.escape(key)}\s*=\s*)(\d+)", re.MULTILINE)
    if pattern.search(section):
        return pattern.sub(lambda match: f"{match.group(1)}{value}", section, count=1)

    trailing_newline = "" if section.endswith("\n") else "\n"
    return f"{section}{trailing_newline}{key} = {value}\n"


def update_section_text(
    content: str,
    species: str,
    total_samples: int,
    total_cells: int,
) -> str:
    """Update numeric fields for a species in TOML content without disturbing comments."""

    section_pattern = re.compile(
        rf"(\[{re.escape(species)}\][\s\S]*?)(?=\n\[|$)",
    )
    match = section_pattern.search(content)
    if not match:
        raise RuntimeError(f"Species '{species}' section not found in species_display.toml")

    section = match.group(1)
    updated_section = replace_or_append(section, "single_cell_datasets", total_samples)
    updated_section = replace_or_append(updated_section, "cells", total_cells)

    start, end = match.span(1)
    return f"{content[:start]}{updated_section}{content[end:]}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize single-cell statistics and update species_display.toml"
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Species key (e.g. zebrafish).",)
    parser.add_argument(
        "--data-root",
        help="Root directory that contains species data folders (default: dash/data)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only compute statistics; do not modify species_display.toml",
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Create a .bak copy of species_display.toml before writing (disabled by default).",
    )
    parser.add_argument(
        "--hugo-toml",
        help="Optional Hugo species_display.toml path to mirror after updating dash/data.",
    )
    parser.add_argument(
        "--no-hugo-mirror",
        action="store_true",
        help="Do not mirror the updated totals into Hugo/data/species_display.toml.",
    )
    return parser.parse_args()


def find_metadata_files(samples_dir: Path) -> Iterable[Tuple[Path, Path]]:
    """Yield pairs of (sample_dir, metadata_file) for valid samples."""

    for study_dir in sorted([d for d in samples_dir.iterdir() if d.is_dir()]):
        if "_Unpub_" in study_dir.name:
            print(f"⏭️  Skip unpublished study: {study_dir}")
            continue
        for sample_dir in sorted([d for d in study_dir.iterdir() if d.is_dir()]):
            metadata_file = sample_dir / f"{sample_dir.name}_metadata.json"
            if metadata_file.exists():
                yield sample_dir, metadata_file
            else:
                print(
                    f"⚠️  Missing metadata: {metadata_file} — sample will be skipped",
                )


def load_sample_stats(
    species: str,
    assay: str,
    samples_dir: Path,
) -> List[SampleStat]:
    """Collect per-sample statistics for a given assay directory."""

    stats: List[SampleStat] = []
    for sample_dir, metadata_file in find_metadata_files(samples_dir):
        try:
            with open(metadata_file, "r", encoding="utf-8") as fh:
                meta = json.load(fh)
        except json.JSONDecodeError as exc:
            print(f"❌ Failed to parse {metadata_file}: {exc}")
            continue

        n_cells = meta.get("n_cells")
        if n_cells is None:
            print(
                f"⚠️  Metadata missing n_cells: {metadata_file} — treated as 0",
            )
            n_cells = 0

        stats.append(
            SampleStat(
                species=species,
                assay=assay,
                study=sample_dir.parent.name,
                sample=sample_dir.name,
                cells=int(n_cells),
            )
        )
    return stats


def write_csv(output_path: Path, rows: Iterable[SampleStat]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["species", "assay", "study", "sample", "cells"])
        for row in rows:
            writer.writerow([row.species, row.assay, row.study, row.sample, row.cells])


def summarize(stats: Iterable[SampleStat]) -> Dict[str, Dict[str, int]]:
    """Return aggregated counts keyed by assay."""

    aggregates: Dict[str, Dict[str, int]] = defaultdict(lambda: {"samples": 0, "cells": 0})
    for stat in stats:
        aggregates[stat.assay]["samples"] += 1
        aggregates[stat.assay]["cells"] += stat.cells
    return aggregates


def print_summary(aggregates: Dict[str, Dict[str, int]]) -> None:
    total_samples = sum(v["samples"] for v in aggregates.values())
    total_cells = sum(v["cells"] for v in aggregates.values())

    print("\n=== Single-cell dataset summary ===")
    for assay, values in sorted(aggregates.items()):
        print(
            f"{assay:>6}: {values['samples']} samples, {values['cells']} cells",
        )
    print("----------------------------------")
    print(f" total: {total_samples} samples, {total_cells} cells\n")


def update_species_display(
    toml_path: Path,
    species: str,
    total_samples: int,
    total_cells: int,
    dry_run: bool,
    backup: bool,
) -> None:
    """Update dash/data/species_display.toml with the provided totals."""

    if not toml_path.exists():
        print(
            f"⚠️  species_display.toml not found: {toml_path}\n"
            "Please copy the latest Hugo/data/species_display.toml into dash/data/ and rerun.",
        )
        return

    original_text: Optional[str] = None
    if tomlkit is not None:
        toml_data = load_toml_file(toml_path)
    else:
        with open(toml_path, "r", encoding="utf-8") as fh:
            original_text = fh.read()
        toml_data = toml.loads(original_text)

    if species not in toml_data:
        print(f"⚠️  Species '{species}' not found in {toml_path}")
        return

    current = toml_data[species]
    old_samples_value = current.get("single_cell_datasets")
    old_cells_value = current.get("cells")
    old_samples = int(old_samples_value) if old_samples_value is not None else 0
    old_cells = int(old_cells_value) if old_cells_value is not None else 0

    print(
        f"Updating species_display.toml for '{species}':"
        f" samples {old_samples} -> {total_samples},"
        f" cells {old_cells} -> {total_cells}"
    )

    if dry_run:
        print("Dry run: no changes written to species_display.toml")
        return

    if backup:
        backup_path = toml_path.with_suffix(".toml.bak")
        shutil.copy2(toml_path, backup_path)
        print(f"Backup written to {backup_path}")

    if tomlkit is not None:
        current["single_cell_datasets"] = make_toml_int(total_samples)
        current["cells"] = make_toml_int(total_cells)
        dump_toml_file(toml_data, toml_path)
    else:
        if original_text is None:  # pragma: no cover - defensive
            raise RuntimeError("Internal error: original TOML text unavailable for update")
        updated_text = update_section_text(
            original_text,
            species,
            total_samples,
            total_cells,
        )
        with open(toml_path, "w", encoding="utf-8") as fh:
            fh.write(updated_text)

    print(f"species_display.toml updated successfully: {toml_path}")


def main() -> None:
    args = parse_args()

    script_dir = Path(__file__).resolve().parent
    default_data_root = script_dir.parent / "data"
    repo_root = script_dir.parents[1]

    species = args.species
    data_root = Path(args.data_root).resolve() if args.data_root else default_data_root
    species_dir = data_root / species

    if not species_dir.exists():
        raise SystemExit(f"Species directory not found: {species_dir}")

    toml_path = data_root / "species_display.toml"

    assays = {
        "scRNA": species_dir / "scRNA",
        "scATAC": species_dir / "scATAC",
    }

    all_stats: List[SampleStat] = []
    for assay, assay_dir in assays.items():
        if not assay_dir.exists():
            print(f"🔍 Skipping {assay}: directory not found ({assay_dir})")
            continue

        stats = load_sample_stats(species, assay, assay_dir)
        if not stats:
            print(f"ℹ️  No valid samples found in {assay_dir}")
            continue

        all_stats.extend(stats)

    if not all_stats:
        print("No valid single-cell samples found; aborting.")
        return

    # Write per-sample CSV next to species data directory.
    csv_path = species_dir / "single_cell_stats.csv"
    write_csv(csv_path, all_stats)
    print(f"Per-sample statistics written to {csv_path}")

    aggregates = summarize(all_stats)
    print_summary(aggregates)

    total_samples = sum(v["samples"] for v in aggregates.values())
    total_cells = sum(v["cells"] for v in aggregates.values())

    update_species_display(
        toml_path=toml_path,
        species=species,
        total_samples=total_samples,
        total_cells=total_cells,
        dry_run=args.dry_run,
        backup=args.backup,
    )

    if args.no_hugo_mirror:
        return

    if args.hugo_toml:
        hugo_toml = Path(args.hugo_toml).resolve()
    else:
        hugo_toml = repo_root / "hugo" / "data" / "species_display.toml"

    if hugo_toml == toml_path:
        return

    if not hugo_toml.exists():
        print(f"ℹ️  Skip Hugo mirror: file not found ({hugo_toml})")
        return

    print(f"Mirroring single-cell totals into Hugo config: {hugo_toml}")
    update_species_display(
        toml_path=hugo_toml,
        species=species,
        total_samples=total_samples,
        total_cells=total_cells,
        dry_run=args.dry_run,
        backup=args.backup,
    )


if __name__ == "__main__":
    main()
