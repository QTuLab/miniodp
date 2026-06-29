#!/usr/bin/env python3
"""
Preprocess JBrowse track files.

Features:
1. Convert .gff3 to .gff3.gz + .gff3.gz.tbi (sort -> bgzip -> tabix)
2. Convert .bed to .bb (clean -> bedToBigBed)
3. Skip already usable format files (.bw, .bam, .bb, etc)

Incremental support:
- If target file exists and is newer than source, skip conversion

Required tools:
- bgzip, tabix (for GFF3)
- bedToBigBed (for BED, requires chrom.sizes)

Usage:
    python3 prepare_tracks.py Danio_rerio
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_APP_ROOT = Path("/opt/miniodp/jbrowse2")
DEFAULT_DATA_ROOT = DEFAULT_APP_ROOT / "data"

SKIP_EXTENSIONS = {
    ".bw", ".bigwig",
    ".bb", ".bigbed",
    ".bam", ".cram",
    ".vcf.gz", ".bcf",
    ".fai", ".tbi", ".csi", ".bai", ".crai", ".gzi",
}

SKIP_COMPOSITE_SUFFIXES = (
    ".gff3.gz",
    ".gff.gz",
    ".gtf.gz",
    ".bed.gz",
    ".vcf.gz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Preprocess JBrowse track files (GFF3 -> gz+tbi, BED -> bb)"
    )
    parser.add_argument("species", help="species directory name or absolute path")
    parser.add_argument(
        "--data-root",
        default=str(DEFAULT_DATA_ROOT),
        help=f"species data root directory, default {DEFAULT_DATA_ROOT}",
    )
    parser.add_argument(
        "--chrom-sizes",
        default=None,
        help="Explicitly specify chrom.sizes file, default: <species>/reference/chrom.sizes",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force reconvert all files",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete original GFF3/GFF/GTF files after successful conversion (BED not deleted as BigBed conversion is lossy)",
    )
    return parser.parse_args()


def resolve_species_dir(species: str, data_root: Path) -> Path:
    candidate = Path(species)
    if candidate.is_dir():
        return candidate.resolve()
    combined = data_root / species
    if combined.is_dir():
        return combined.resolve()
    sys.exit(f"Error: Not found:species directory '{species}' (data root: {data_root})")


def is_newer(target: Path, source: Path) -> bool:
    """Check if target exists and is newer than source"""
    if not target.exists():
        return False
    return target.stat().st_mtime >= source.stat().st_mtime


def load_chrom_sizes(path: Path) -> dict[str, int]:
    if not path.exists():
        sys.exit(f"Error: Not found: chrom.sizes ({path}), please run prepare_reference.py first")
    sizes: dict[str, int] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def should_skip_file(path: Path) -> bool:
    """Check if file should be skipped (already usable format or index file)"""
    name_lower = path.name.lower()
    for suffix in SKIP_COMPOSITE_SUFFIXES:
        if name_lower.endswith(suffix):
            return True
    suffix = path.suffix.lower()
    return suffix in SKIP_EXTENSIONS


def sort_gff(gff_path: Path) -> Path:
    """Sort GFF3 file, return temp file path"""
    headers: list[str] = []
    records: list[tuple[str, int, str]] = []
    with gff_path.open("r") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                headers.append(line)
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            try:
                start = int(fields[3])
            except ValueError:
                continue
            records.append((fields[0], start, line))
    records.sort(key=lambda item: (item[0], item[1]))
    fd, tmp_name = tempfile.mkstemp(suffix=".gff3")
    os.close(fd)
    tmp = Path(tmp_name)
    with tmp.open("w") as out:
        out.writelines(headers)
        for _, _, line in records:
            out.write(line)
    return tmp


def convert_gff3(gff_path: Path, force: bool, clean: bool) -> tuple[bool, str]:
    """
    Convert GFF3 file to bgzip compressed format and build index.
    Return (success, message)
    """
    gz_path = gff_path.with_suffix(gff_path.suffix + ".gz")
    tbi_path = Path(str(gz_path) + ".tbi")

    if not force and is_newer(gz_path, gff_path) and is_newer(tbi_path, gff_path):
        return True, f"Skip: {gz_path.name} exists and is newer"

    try:
        sorted_tmp = sort_gff(gff_path)
        with gz_path.open("wb") as gz_out:
            subprocess.run(
                ["bgzip", "-c", str(sorted_tmp)],
                stdout=gz_out,
                check=True,
            )
        subprocess.run(
            ["tabix", "-f", "-p", "gff", str(gz_path)],
            check=True,
        )
        sorted_tmp.unlink(missing_ok=True)
        msg = f"Conversion done: {gff_path.name} -> {gz_path.name}"
        if clean:
            gff_path.unlink()
            msg += " (original file deleted)"
        return True, msg
    except subprocess.CalledProcessError as e:
        return False, f"Conversion failed: {gff_path.name} ({e})"
    except Exception as e:
        return False, f"Conversion failed: {gff_path.name} ({e})"


def sanitize_bed(bed_path: Path, chrom_sizes: dict[str, int]) -> tuple[Path | None, dict[str, int]]:
    """Clean BED file, return (temp file path, statistics)"""
    stats = {"total": 0, "kept": 0, "clipped": 0, "skipped": 0}
    records: list[tuple[str, int, int, str]] = []
    with bed_path.open() as fh:
        for line in fh:
            raw = line.strip()
            if not raw or raw.startswith(("#", "track", "browser")):
                continue
            stats["total"] += 1
            parts = raw.split()
            if len(parts) < 3:
                stats["skipped"] += 1
                continue
            chrom = parts[0]
            if chrom not in chrom_sizes:
                stats["skipped"] += 1
                continue
            try:
                start = int(float(parts[1]))
                end = int(float(parts[2]))
            except ValueError:
                stats["skipped"] += 1
                continue
            if end <= start:
                stats["skipped"] += 1
                continue
            chrom_len = chrom_sizes[chrom]
            clipped = False
            if start < 0:
                start = 0
                clipped = True
            if end > chrom_len:
                end = chrom_len
                clipped = True
            if start >= chrom_len or start >= end:
                stats["skipped"] += 1
                continue
            name = parts[3] if len(parts) >= 4 else "."
            if clipped:
                stats["clipped"] += 1
            records.append((chrom, start, end, name))
    if not records:
        return None, stats
    records.sort(key=lambda item: (item[0], item[1]))
    fd, tmp_name = tempfile.mkstemp(suffix=".bed")
    os.close(fd)
    tmp_path = Path(tmp_name)
    with tmp_path.open("w") as out:
        for chrom, start, end, name in records:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")
            stats["kept"] += 1
    return tmp_path, stats


def convert_bed(bed_path: Path, chrom_sizes_path: Path, chrom_sizes: dict[str, int], force: bool) -> tuple[bool, str]:
    """
    Convert BED file to BigBed format.
    Return (success, message)
    """
    bb_path = bed_path.with_suffix(".bb")

    if not force and is_newer(bb_path, bed_path):
        return True, f"Skip: {bb_path.name} exists and is newer"

    sorted_tmp, stats = sanitize_bed(bed_path, chrom_sizes)
    if sorted_tmp is None:
        return False, f"Skip: {bed_path.name} no valid records"

    try:
        subprocess.run(
            ["bedToBigBed", "-type=bed4", str(sorted_tmp), str(chrom_sizes_path), str(bb_path)],
            check=True,
            capture_output=True,
        )
        msg = f"Conversion done: {bed_path.name} -> {bb_path.name}"
        if stats["clipped"]:
            msg += f" (clipped {stats['clipped']} out-of-bounds records)"
        return True, msg
    except subprocess.CalledProcessError as e:
        return False, f"Conversion failed: {bed_path.name} ({e.stderr.decode() if e.stderr else e})"
    finally:
        if sorted_tmp:
            sorted_tmp.unlink(missing_ok=True)


def find_track_files(species_dir: Path) -> list[Path]:
    """Find files to process under tracks/ directory"""
    tracks_dir = species_dir / "tracks"
    if not tracks_dir.is_dir():
        return []

    files: list[Path] = []
    for path in tracks_dir.rglob("*"):
        if not path.is_file():
            continue
        if path.name.startswith("."):
            continue
        if should_skip_file(path):
            continue
        suffix = path.suffix.lower()
        if suffix in (".gff3", ".gff", ".gtf", ".bed"):
            files.append(path)
    return sorted(files)


def main() -> None:
    args = parse_args()
    data_root = Path(args.data_root).resolve()
    species_dir = resolve_species_dir(args.species, data_root)
    force = args.force
    clean = args.clean

    chrom_sizes_path = Path(args.chrom_sizes) if args.chrom_sizes else species_dir / "reference" / "chrom.sizes"
    chrom_sizes: dict[str, int] | None = None

    track_files = find_track_files(species_dir)
    if not track_files:
        print("Note: no track files found for conversion (.gff3, .gff, .gtf, .bed)")
        return

    print(f"Found {len(track_files)} files to process")
    summary = {"success": 0, "skipped": 0, "failed": 0, "cleaned": 0}

    for track_file in track_files:
        suffix = track_file.suffix.lower()

        if suffix in (".gff3", ".gff", ".gtf"):
            success, msg = convert_gff3(track_file, force, clean)
            if clean and "original file deleted" in msg:
                summary["cleaned"] += 1
        elif suffix == ".bed":
            if chrom_sizes is None:
                chrom_sizes = load_chrom_sizes(chrom_sizes_path)
            success, msg = convert_bed(track_file, chrom_sizes_path, chrom_sizes, force)
        else:
            continue

        print(f"  {msg}")
        if "skipped" in msg:
            summary["skipped"] += 1
        elif success:
            summary["success"] += 1
        else:
            summary["failed"] += 1

    result_msg = f"\nDone: converted {summary['success']}, skipped {summary['skipped']}, failed {summary['failed']}"
    if summary["cleaned"]:
        result_msg += f", cleaned {summary['cleaned']} original files"
    print(result_msg)
    if summary["failed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
