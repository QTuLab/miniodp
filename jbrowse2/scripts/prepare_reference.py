#!/usr/bin/env python3
"""
Prepare JBrowse reference genome for a species.

Features:
1. Generate index (.fai) and chrom.sizes for reference/*.fa
2. Call jbrowse add-assembly to register reference genome

Incremental support:
- If .fai and chrom.sizes already existsand is newer than source, skip generation
- Skip registration if assembly already in config.json

Usage:
    python3 prepare_reference.py Danio_rerio
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

DEFAULT_APP_ROOT = Path("/opt/miniodp/jbrowse2")
DEFAULT_DATA_ROOT = DEFAULT_APP_ROOT / "data"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare JBrowse reference genome for a species"
    )
    parser.add_argument("species", help="species directory name or absolute path")
    parser.add_argument(
        "--assembly-name",
        default=None,
        help="Register to JBrowse assembly name (default: species directory name)",
    )
    parser.add_argument(
        "--app-root",
        default=str(DEFAULT_APP_ROOT),
        help=f"JBrowse app directory, default {DEFAULT_APP_ROOT}",
    )
    parser.add_argument(
        "--data-root",
        default=str(DEFAULT_DATA_ROOT),
        help=f"species data root directory, default {DEFAULT_DATA_ROOT}",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force regenerate index and register assembly",
    )
    return parser.parse_args()


def run_cmd(cmd: list[str], *, cwd: Path | None = None) -> None:
    print(f"$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=cwd)


def path_for_jbrowse(target: Path, app_root: Path) -> str:
    try:
        return str(target.resolve().relative_to(app_root))
    except ValueError:
        return str(target)


def resolve_species_dir(species: str, data_root: Path) -> Path:
    candidate = Path(species)
    if candidate.is_dir():
        return candidate.resolve()
    combined = data_root / species
    if combined.is_dir():
        return combined.resolve()
    sys.exit(f"Error: Not found:species directory '{species}' (data root: {data_root})")


def pick_single_file(directory: Path, patterns: tuple[str, ...]) -> Path:
    files: list[Path] = []
    for pattern in patterns:
        files.extend(directory.glob(pattern))
    files = [f for f in files if f.is_file()]
    if not files:
        sys.exit(f"Error: {directory}  cannot find matching {patterns} files")
    if len(files) > 1:
        pretty = "\n  - ".join(str(f) for f in files)
        sys.exit(f"Error: {directory} has multiple candidate files, please keep only one:\n  - {pretty}")
    return files[0]


def is_newer(target: Path, source: Path) -> bool:
    """Check if target exists and is newer than source"""
    if not target.exists():
        return False
    return target.stat().st_mtime >= source.stat().st_mtime


def load_existing_assemblies(app_root: Path) -> set[str]:
    """Read registered assembly names from config.json"""
    config_path = app_root / "config.json"
    if not config_path.exists():
        return set()
    try:
        data = json.loads(config_path.read_text())
        return {asm.get("name") for asm in data.get("assemblies", []) if asm.get("name")}
    except (json.JSONDecodeError, KeyError):
        return set()


def ensure_fai(fasta: Path, force: bool) -> tuple[Path, bool]:
    """Generate FASTA index, return (fai path, whether newly generated)"""
    fai = fasta.with_suffix(fasta.suffix + ".fai")
    if not force and is_newer(fai, fasta):
        print(f"Skip: {fai.name} exists and is newer")
        return fai, False
    print(f"Generating FASTA index -> {fai}")
    run_cmd(["samtools", "faidx", str(fasta)])
    if not fai.exists():
        sys.exit(f"Error: samtools faidx did not generate {fai}")
    return fai, True


def write_chrom_sizes(fai: Path, output: Path, force: bool) -> bool:
    """Generate chrom.sizes, return whether newly generated"""
    if not force and is_newer(output, fai):
        print(f"Skip: {output.name} exists and is newer")
        return False
    print(f"Generating chrom.sizes -> {output}")
    with fai.open("r") as src, output.open("w") as dest:
        for line in src:
            if not line.strip():
                continue
            fields = line.split("\t")
            dest.write(f"{fields[0]}\t{fields[1]}\n")
    return True


def main() -> None:
    args = parse_args()
    app_root = Path(args.app_root).resolve()
    data_root = Path(args.data_root).resolve()
    species_dir = resolve_species_dir(args.species, data_root)
    force = args.force

    reference_dir = species_dir / "reference"
    if not reference_dir.is_dir():
        sys.exit(f"Error: Missing reference directory {reference_dir}")

    fasta_path = pick_single_file(reference_dir, ("*.fa", "*.fasta"))
    fai_path, fai_created = ensure_fai(fasta_path, force)
    chrom_sizes = reference_dir / "chrom.sizes"
    chrom_created = write_chrom_sizes(fai_path, chrom_sizes, force)
    refname_aliases = reference_dir / "refname_aliases.tsv"

    assembly = args.assembly_name or species_dir.name
    existing_assemblies = load_existing_assemblies(app_root)

    if not force and assembly in existing_assemblies:
        print(f"Skip: assembly '{assembly}' already registered in config.json")
        print(f"Done: species={species_dir.name}, assembly={assembly} (no update needed)")
        return

    print(f"Registering assembly: {assembly}")
    jbrowse_fasta = path_for_jbrowse(fasta_path, app_root)
    add_assembly_cmd = [
        "jbrowse",
        "add-assembly",
        jbrowse_fasta,
        "--name",
        assembly,
        "--load",
        "inPlace",
    ]
    if refname_aliases.exists():
        add_assembly_cmd.extend(
            [
                "--refNameAliases",
                path_for_jbrowse(refname_aliases, app_root),
            ]
        )

    if assembly in existing_assemblies:
        print(f"Note: assembly '{assembly}' already exists, will use --overwrite to override")
        run_cmd(add_assembly_cmd + ["--overwrite"], cwd=app_root)
    else:
        run_cmd(add_assembly_cmd, cwd=app_root)

    print(f"Done: species={species_dir.name}, assembly={assembly}")


if __name__ == "__main__":
    main()
