#!/usr/bin/env python3
"""
Read MetaData.csv and generate JBrowse add-track script.

Features:
1. Parse track metadata from MetaData.csv
2. Auto-match data files under tracks/ directory
3. Generate batch execution bash script

Incremental support:
- Check existing tracks in config.json and skip automatically

File matching priority:
- Prioritize converted formats (.gff3.gz > .gff3, .bb > .bed)
- If only raw format found, prompt user to run prepare_tracks.py first

Usage:
    python3 generate_loader.py Danio_rerio
    bash data/Danio_rerio/run_add_tracks.sh
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import sys
from collections import defaultdict
from pathlib import Path

DEFAULT_APP_ROOT = "/opt/miniodp/jbrowse2"
DEFAULT_DATA_ROOT = os.path.join(DEFAULT_APP_ROOT, "data")
TRACK_DIR_NAME = "tracks"

SKIP_INDEX_SUFFIXES = {".bai", ".crai", ".csi", ".tbi", ".fai", ".gzi"}
SKIP_FILE_NAMES = {"MetaData.csv", "MetaData2.csv"}

RAW_FORMATS = {".gff3", ".gff", ".gtf", ".bed"}
DIRECT_LOAD_RAW_FORMATS = {".gtf"}

CONVERTED_FORMATS = {
    ".gff3.gz": ".gff3",
    ".gff.gz": ".gff",
    ".bb": ".bed",
}

UNSUPPORTED_AUTOLOAD_FORMATS = {
    ".gtf.gz": ".gtf",
}

COMPOSITE_SUFFIXES = (
    ".bam.bw",
    ".bam.bai",
    ".cram.crai",
    ".bed.gz",
    ".gff3.gz",
    ".gff.gz",
    ".gtf.gz",
    ".vcf.gz",
    ".fa.gz",
    ".fasta.gz",
)

OPTIONAL_NAME_SUFFIXES = (
    ".final",
    ".merge",
    ".merged",
    ".dedup",
    ".filtered",
    ".sorted",
    "_final",
    "_dedup",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read MetaData.csv and generate JBrowse add-track script"
    )
    parser.add_argument("species", help="species directory name or absolute path")
    parser.add_argument(
        "--data-root",
        default=DEFAULT_DATA_ROOT,
        help=f"species data root directory (default: {DEFAULT_DATA_ROOT})",
    )
    parser.add_argument(
        "--app-root",
        default=DEFAULT_APP_ROOT,
        help=f"jbrowse app directory (default: {DEFAULT_APP_ROOT})",
    )
    parser.add_argument(
        "--assembly-name",
        default=None,
        help="assembly name used in jbrowse, defaults to species directory name",
    )
    parser.add_argument(
        "--metadata-name",
        default="MetaData.csv",
        help="MetaData filename (default: MetaData.csv)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Exported bash file path, defaults to <species>/run_add_tracks.sh",
    )
    return parser.parse_args()


def resolve_species_dir(species: str, data_root: str) -> Path:
    candidate = Path(species)
    if candidate.is_dir():
        return candidate.resolve()
    base = Path(data_root) / species
    if base.is_dir():
        return base.resolve()
    sys.exit(f"Error: Not found:species directory '{species}' (data root: {data_root})")


def _iter_metadata_lines(fh):
    for line in fh:
        stripped = line.lstrip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue
        yield line


def read_metadata(metadata_path: Path) -> list[dict[str, str]]:
    if not metadata_path.exists():
        sys.exit(f"Error: Metadata file not found {metadata_path}")
    rows: list[dict[str, str]] = []
    with metadata_path.open("r", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(_iter_metadata_lines(fh))
        missing = {"label"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"Error: MetaData.csv missing fields: {', '.join(sorted(missing))}")
        for row in reader:
            if not row.get("label"):
                continue
            rows.append({k: v.strip() if isinstance(v, str) else v for k, v in row.items()})
    if not rows:
        sys.exit("Error: MetaData.csv has no valid records")
    return rows


def path_for_jbrowse(target: Path, app_root: Path) -> str:
    try:
        return str(target.resolve().relative_to(app_root))
    except ValueError:
        return str(target)


def get_file_format(path: Path) -> str:
    """Get file format (supports compound suffix)"""
    name_lower = path.name.lower()
    for suffix in COMPOSITE_SUFFIXES:
        if name_lower.endswith(suffix):
            return suffix
    return path.suffix.lower()


def strip_known_extensions(name: str) -> str:
    lower = name.lower()
    for suffix in COMPOSITE_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    if "." in name:
        return name.rsplit(".", 1)[0]
    return name


def label_keys_from_name(name: str) -> list[str]:
    base = strip_known_extensions(name)
    keys: list[str] = []
    seen: set[str] = set()

    def add(value: str) -> None:
        trimmed = value.strip()
        if not trimmed or trimmed in seen:
            return
        seen.add(trimmed)
        keys.append(trimmed)

    add(base)
    queue = [base]
    while queue:
        current = queue.pop()
        lowered = current.lower()
        for suffix in OPTIONAL_NAME_SUFFIXES:
            if lowered.endswith(suffix):
                trimmed = current[: -len(suffix)]
                if trimmed.strip() and trimmed not in seen:
                    add(trimmed)
                    queue.append(trimmed)
    return keys


def should_skip_file(path: Path) -> bool:
    if path.name in SKIP_FILE_NAMES:
        return True
    if path.name.startswith("."):
        return True
    suffixes = [s.lower() for s in path.suffixes]
    if suffixes and suffixes[-1] in SKIP_INDEX_SUFFIXES:
        return True
    return False


def build_label_index(species_dir: Path) -> dict[str, list[Path]]:
    roots: list[Path] = []
    tracks_dir = species_dir / TRACK_DIR_NAME
    if tracks_dir.is_dir():
        roots.append(tracks_dir)
    else:
        roots.append(species_dir)

    index: dict[str, list[Path]] = defaultdict(list)
    for root in roots:
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if should_skip_file(path):
                continue
            for key in label_keys_from_name(path.name):
                index[key].append(path)
    return index


def select_track_path(
    label: str, index: dict[str, list[Path]], species_dir: Path
) -> tuple[Path, bool]:
    """
    Select track file path.
    Returns (file path, whether raw format needs conversion)
    """
    candidates = index.get(label)
    if not candidates:
        sys.exit(
            f"Error: not found label '{label}' corresponding data file; please confirm filename = label"
        )

    unique_paths = list({p.resolve() for p in candidates})

    converted_files: list[Path] = []
    direct_load_raw_files: list[Path] = []
    raw_files: list[Path] = []
    unsupported_autoload_files: list[Path] = []
    other_files: list[Path] = []

    for path in unique_paths:
        fmt = get_file_format(path)
        if fmt in CONVERTED_FORMATS:
            converted_files.append(path)
        elif fmt in DIRECT_LOAD_RAW_FORMATS:
            direct_load_raw_files.append(path)
        elif fmt in RAW_FORMATS:
            raw_files.append(path)
        elif fmt in UNSUPPORTED_AUTOLOAD_FORMATS:
            unsupported_autoload_files.append(path)
        else:
            other_files.append(path)

    if converted_files:
        if len(converted_files) > 1:
            formatted = [str(p.relative_to(species_dir)) if p.is_relative_to(species_dir) else str(p) for p in converted_files]
            sys.exit(f"Error: label '{label}' matched multiple converted files: {', '.join(formatted)}")
        return converted_files[0], False

    if direct_load_raw_files:
        if len(direct_load_raw_files) > 1:
            formatted = [str(p.relative_to(species_dir)) if p.is_relative_to(species_dir) else str(p) for p in direct_load_raw_files]
            sys.exit(f"Error: label '{label}' matched multiple raw files: {', '.join(formatted)}")
        return direct_load_raw_files[0], False

    if unsupported_autoload_files:
        formatted = [str(p.relative_to(species_dir)) if p.is_relative_to(species_dir) else str(p) for p in unsupported_autoload_files]
        sys.exit(
            f"Error: label '{label}' only matched unsupported autoload format: {', '.join(formatted)}; "
            "please keep the raw .gtf file or convert the annotation to .gff3/.gff3.gz"
        )

    if other_files:
        if len(other_files) > 1:
            formatted = [str(p.relative_to(species_dir)) if p.is_relative_to(species_dir) else str(p) for p in other_files]
            sys.exit(f"Error: label '{label}' matched multiple files: {', '.join(formatted)}")
        return other_files[0], False

    if raw_files:
        if len(raw_files) > 1:
            formatted = [str(p.relative_to(species_dir)) if p.is_relative_to(species_dir) else str(p) for p in raw_files]
            sys.exit(f"Error: label '{label}' matched multiple raw files: {', '.join(formatted)}")
        return raw_files[0], True

    sys.exit(f"Error: not found label '{label}' corresponding available data file")


def build_command(
    row: dict[str, str], file_path: Path, assembly: str, app_root: Path
) -> tuple[str, str]:
    label = row.get("label", "").strip()

    if not file_path.exists():
        sys.exit(f"Error: Not found:data file {file_path} (label={label})")

    possible_labels = set(label_keys_from_name(file_path.name))
    if label not in possible_labels:
        print(
            f"Note: label '{label}' does not exactly match filename, will import by label",
            file=sys.stderr,
        )

    category = row.get("category", "").strip() or file_path.parent.name
    description = row.get("description", "").strip()

    metadata_payload = {
        k: v for k, v in row.items() if k not in {"relative_path"} and v
    }
    metadata_json = json.dumps({"metadata": metadata_payload}, ensure_ascii=False)
    metadata_arg = shlex.quote(metadata_json)

    jbrowse_path = path_for_jbrowse(file_path, app_root)

    cmd_lines = [
        f"if track_exists {shlex.quote(label)}; then",
        f"  echo \"Skip existing track: {label}\"",
        "else",
        f"  jbrowse add-track {shlex.quote(jbrowse_path)} \\",
        f"    --name {shlex.quote(label)} \\",
        f"    --assemblyNames {shlex.quote(assembly)} \\",
        f"    --category {shlex.quote(category)} \\",
        f"    --load inPlace \\",
    ]
    if description:
        cmd_lines.append(f"    --description {shlex.quote(description)} \\")
    cmd_lines.append(f"    --config {metadata_arg}")
    cmd_lines.append("fi")
    return label, "\n".join(cmd_lines)


def main() -> None:
    args = parse_args()
    species_dir = resolve_species_dir(args.species, args.data_root)
    assembly = args.assembly_name or species_dir.name
    metadata_path = species_dir / args.metadata_name
    rows = read_metadata(metadata_path)
    app_root_resolved = Path(args.app_root).resolve()
    if not app_root_resolved.exists():
        sys.exit(f"Error: Not found: JBrowse directory {app_root_resolved}")
    track_index = build_label_index(species_dir)

    seen_labels: set[str] = set()
    commands: list[str] = []
    needs_conversion: list[str] = []

    for row in rows:
        label = (row.get("label") or "").strip()
        if label in seen_labels:
            sys.exit(f"Error: Duplicate label in MetaData.csv: {label}")
        seen_labels.add(label)
        file_path, is_raw = select_track_path(label, track_index, species_dir)
        if is_raw:
            needs_conversion.append(f"  - {label}: {file_path.name}")
            continue
        _, cmd = build_command(row, file_path, assembly, app_root_resolved)
        commands.append(cmd)

    if needs_conversion:
        print(
            f"\nWarning: following {len(needs_conversion)} tracks need format conversion first, please run prepare_tracks.py:",
            file=sys.stderr,
        )
        for item in needs_conversion:
            print(item, file=sys.stderr)
        print(f"\n  python3 jbrowse2/scripts/prepare_tracks.py {species_dir.name}\n", file=sys.stderr)

    if not commands:
        if needs_conversion:
            print("No script generated: all new tracks need format conversion first", file=sys.stderr)
        else:
            print("No script generated", file=sys.stderr)
        return

    output_path = Path(args.output) if args.output else species_dir / "run_add_tracks.sh"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    script_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'SCRIPT_PATH="${BASH_SOURCE[0]}"',
        'while [ -L "$SCRIPT_PATH" ]; do',
        '  SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"',
        '  SCRIPT_PATH="$(readlink "$SCRIPT_PATH")"',
        '  [[ "$SCRIPT_PATH" != /* ]] && SCRIPT_PATH="$SCRIPT_DIR/$SCRIPT_PATH"',
        'done',
        'SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"',
        'APP_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/../.." && pwd)"',
        'APP_ROOT="${JBROWSE_APP_ROOT:-$APP_ROOT_DEFAULT}"',
        'if [ ! -d "$APP_ROOT" ]; then',
        '  echo "Error: Not found: JBrowse directory $APP_ROOT" >&2',
        '  exit 1',
        'fi',
        'cd "$APP_ROOT"',
        '',
        'track_exists() {',
        '  local track_name="$1"',
        '  python3 - "$APP_ROOT/config.json" "$track_name" <<'"'"'PY'"'"'',
        'import json',
        'import sys',
        'from pathlib import Path',
        '',
        'config_path = Path(sys.argv[1])',
        'track_name = sys.argv[2]',
        'if not config_path.exists():',
        '    raise SystemExit(1)',
        '',
        'try:',
        '    config = json.loads(config_path.read_text())',
        'except Exception as exc:',
        '    print(f"Error: failed to read {config_path}: {exc}", file=sys.stderr)',
        '    raise SystemExit(2)',
        '',
        'for track in config.get("tracks", []):',
        '    if str(track.get("name", "")) == track_name:',
        '        raise SystemExit(0)',
        '',
        'raise SystemExit(1)',
        'PY',
        '  local status=$?',
        '  if [ "$status" -eq 2 ]; then',
        '    echo "Error: config.json is not readable; abort add-track to avoid false duplicate checks" >&2',
        '    exit 1',
        '  fi',
        '  return "$status"',
        '}',
        '',
        f'echo "StartImporting {len(commands)} tracks (assembly={assembly})"',
        '',
    ]

    for cmd in commands:
        script_lines.append(cmd)
        script_lines.append("")
    script_lines.append('echo "Done with add-track, remember to run jbrowse text-index if needed"')

    output_path.write_text("\n".join(script_lines), encoding="utf-8")
    output_path.chmod(0o755)

    print(
        f"Generated: {output_path} (added {len(commands)} records)"
    )
    if needs_conversion:
        print(f"Note: still have {len(needs_conversion)} tracks pending conversion before running this script again")


if __name__ == "__main__":
    main()
