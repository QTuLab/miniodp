#!/usr/bin/env python3
"""Reorder JBrowse tracks in config.json to match MetaData.csv label order."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to JBrowse config.json")
    parser.add_argument("--metadata", required=True, help="Path to MetaData.csv")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def read_metadata_order(path: Path) -> dict[str, int]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return {row["label"]: index for index, row in enumerate(rows)}


def main() -> int:
    args = parse_args()
    config_path = Path(args.config)
    metadata_path = Path(args.metadata)
    config = json.loads(config_path.read_text())
    tracks = config.get("tracks", [])
    order = read_metadata_order(metadata_path)

    annotated = []
    unknown = []
    for index, track in enumerate(tracks):
        label = str(track.get("name", ""))
        if label in order:
            annotated.append((order[label], track))
        else:
            unknown.append((index, track))

    annotated.sort(key=lambda item: item[0])
    reordered = [track for _, track in annotated] + [track for _, track in unknown]
    changed = reordered != tracks
    print(
        f"tracks={len(tracks)} matched={len(annotated)} unmatched={len(unknown)} changed={changed}"
    )

    if args.dry_run or not changed:
        return 0

    config["tracks"] = reordered
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
