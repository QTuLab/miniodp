#!/usr/bin/env python3
"""Materialize and build SequenceServer BLAST databases."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path


DOCKER_IMAGE = "wurmlab/sequenceserver:latest"
INDEX_SUFFIXES = (
    ".ndb",
    ".nhr",
    ".nin",
    ".njs",
    ".nog",
    ".nos",
    ".not",
    ".nsq",
    ".ntf",
    ".nto",
    ".pdb",
    ".phr",
    ".pin",
    ".pjs",
    ".pog",
    ".pos",
    ".pot",
    ".psq",
    ".ptf",
    ".pto",
)


@dataclass(frozen=True)
class BlastRecord:
    species_key: str
    data_role: str
    file_name: str
    database_title: str
    taxid: int
    dbtype: str
    source_kind: str
    source: str
    parse_seqids: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    common_parser = argparse.ArgumentParser(add_help=False)

    common_parser.add_argument(
        "--config",
        type=Path,
        default=project_dir / "config" / "blastdb_manifest.csv",
        help="BLAST database manifest CSV.",
    )
    common_parser.add_argument(
        "--data-dir",
        type=Path,
        default=script_dir.parent / "data",
        help="Target data directory for FASTA files and BLAST indexes.",
    )
    common_parser.add_argument(
        "--species",
        action="append",
        default=[],
        help="Species key to include. Repeat to select multiple species.",
    )
    common_parser.add_argument(
        "--data-role",
        action="append",
        default=[],
        help="Data role to include, such as genome, cdna, or protein.",
    )
    common_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing FASTA files during materialize.",
    )
    common_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned actions without modifying files.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("list", help="List manifest entries.", parents=[common_parser])
    subparsers.add_parser(
        "materialize",
        help="Download or copy managed FASTA files into the data directory.",
        parents=[common_parser],
    )
    subparsers.add_parser(
        "build",
        help="Run makeblastdb for FASTA files present in the data directory.",
        parents=[common_parser],
    )
    return parser.parse_args()


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def load_records(config_path: Path) -> list[BlastRecord]:
    with config_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        records = [
            BlastRecord(
                species_key=row["species_key"].strip(),
                data_role=row["data_role"].strip(),
                file_name=row["file_name"].strip(),
                database_title=row["database_title"].strip(),
                taxid=int(row["taxid"]),
                dbtype=row["dbtype"].strip(),
                source_kind=row["source_kind"].strip(),
                source=row["source"].strip(),
                parse_seqids=parse_bool(row["parse_seqids"]),
            )
            for row in reader
        ]
    if not records:
        raise ValueError(f"No BLAST records found in {config_path}")
    return records


def filter_records(records: list[BlastRecord], args: argparse.Namespace) -> list[BlastRecord]:
    species_filter = set(args.species)
    role_filter = set(args.data_role)

    filtered = []
    for record in records:
        if species_filter and record.species_key not in species_filter:
            continue
        if role_filter and record.data_role not in role_filter:
            continue
        filtered.append(record)
    return filtered


def repo_root_from_config(config_path: Path) -> Path:
    return config_path.resolve().parents[2]


def resolve_local_source(repo_root: Path, source: str) -> Path:
    source_path = Path(source)
    if source_path.is_absolute():
        return source_path
    return repo_root / source_path


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def materialize_url(url: str, target_path: Path, *, force: bool, dry_run: bool) -> None:
    if target_path.exists() and not force:
        print(f"Skip existing: {target_path}")
        return
    print(f"Fetch: {url} -> {target_path}")
    if dry_run:
        return

    ensure_parent(target_path)
    with tempfile.NamedTemporaryFile(delete=False, dir=target_path.parent) as tmp_handle:
        tmp_path = Path(tmp_handle.name)

    try:
        with urllib.request.urlopen(url) as response:
            if url.endswith(".gz"):
                with gzip.GzipFile(fileobj=response) as gz_handle, tmp_path.open("wb") as out_handle:
                    shutil.copyfileobj(gz_handle, out_handle)
            else:
                with tmp_path.open("wb") as out_handle:
                    shutil.copyfileobj(response, out_handle)
        tmp_path.replace(target_path)
        target_path.chmod(0o644)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def materialize_local(source_path: Path, target_path: Path, *, force: bool, dry_run: bool) -> None:
    if not source_path.exists():
        print(f"Skip missing local source: {source_path}")
        return
    if target_path.exists() and not force:
        print(f"Skip existing: {target_path}")
        return
    print(f"Copy: {source_path} -> {target_path}")
    if dry_run:
        return
    ensure_parent(target_path)
    shutil.copyfile(source_path, target_path)
    target_path.chmod(0o644)


def materialize_records(
    records: list[BlastRecord],
    *,
    config_path: Path,
    data_dir: Path,
    force: bool,
    dry_run: bool,
) -> None:
    repo_root = repo_root_from_config(config_path)
    for record in records:
        if not record.source_kind or not record.source:
            print(f"Skip unmanaged source: {record.file_name}")
            continue

        target_path = data_dir / record.file_name
        if record.source_kind == "url":
            materialize_url(record.source, target_path, force=force, dry_run=dry_run)
            continue
        if record.source_kind == "local":
            source_path = resolve_local_source(repo_root, record.source)
            materialize_local(source_path, target_path, force=force, dry_run=dry_run)
            continue
        raise ValueError(f"Unsupported source_kind for {record.file_name}: {record.source_kind}")


def remove_existing_indexes(fasta_path: Path, *, dry_run: bool) -> None:
    for suffix in INDEX_SUFFIXES:
        index_path = fasta_path.with_name(fasta_path.name + suffix)
        if not index_path.exists():
            continue
        print(f"Remove stale index: {index_path}")
        if not dry_run:
            index_path.unlink()


def build_records(records: list[BlastRecord], *, data_dir: Path, dry_run: bool) -> None:
    uid_gid = f"{os.getuid()}:{os.getgid()}"
    for record in records:
        fasta_path = data_dir / record.file_name
        if not fasta_path.exists():
            print(f"Skip missing FASTA: {fasta_path}")
            continue

        remove_existing_indexes(fasta_path, dry_run=dry_run)
        cmd = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{data_dir.resolve()}:/db",
            "-u",
            uid_gid,
            "--entrypoint",
            "makeblastdb",
            DOCKER_IMAGE,
            "-in",
            f"/db/{record.file_name}",
            "-dbtype",
            record.dbtype,
            "-title",
            record.database_title,
            "-taxid",
            str(record.taxid),
        ]
        if record.parse_seqids:
            cmd.append("-parse_seqids")
        cmd.extend(
            [
                "-out",
                f"/db/{record.file_name}",
            ]
        )
        print("$ " + " ".join(cmd))
        if not dry_run:
            subprocess.run(cmd, check=True)


def list_records(records: list[BlastRecord]) -> None:
    for record in records:
        source_label = record.source_kind or "manual"
        print(
            "\t".join(
                [
                    record.species_key,
                    record.data_role,
                    record.file_name,
                    record.dbtype,
                    source_label,
                    record.database_title,
                ]
            )
        )


def main() -> int:
    args = parse_args()
    records = filter_records(load_records(args.config.resolve()), args)
    if not records:
        print("No BLAST records matched the requested filters.", file=sys.stderr)
        return 1

    data_dir = args.data_dir.resolve()
    if args.command == "list":
        list_records(records)
        return 0
    if args.command == "materialize":
        materialize_records(
            records,
            config_path=args.config.resolve(),
            data_dir=data_dir,
            force=args.force,
            dry_run=args.dry_run,
        )
        return 0
    if args.command == "build":
        build_records(records, data_dir=data_dir, dry_run=args.dry_run)
        return 0
    raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
