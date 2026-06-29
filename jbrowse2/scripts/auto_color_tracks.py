#!/usr/bin/env python3
"""
Auto-assign colors to JBrowse 2 tracks based on category metadata.

Usage:
    python jbrowse2/scripts/auto_color_tracks.py --config /opt/miniodp/jbrowse2/config.json

The script inspects each track's metadata/category and applies the unified palette:
    - RNA    (#FC8D62)  -> categories containing "rna"
    - ATAC   (#4EB3D3)  -> categories containing "atac"
    - ChIP   (#E78AC3)  -> categories containing "chip" or "histone"
    - Other  (#74C476)  -> categories containing "other" (fallback)

Annotation tracks keep the default styling to preserve CDS/UTR colors.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

# # miniodp v1
# DEFAULT_COLOR_RULES = [
#     ("rna-seq", "#FC8D62"),
#     ("atac-seq", "#4EB3D3"),
#     ("chip-seq", "#E78AC3"),
#     ("derived", "#FFD92F"),
#     ("otheromics", "#74C476"),
# ]

# # NPG (Nature Publishing Group) Inspired
# DEFAULT_COLOR_RULES = [
#     ("rna-seq", "#E64B35"),  # Brick red (emphasis)
#     ("atac-seq", "#4DBBD5"), # Lake blue (clear structure)
#     ("chip-seq", "#00A087"), # Blue-green (stable)
#     ("derived", "#3C5488"),  # Dark blue (good for background/auxiliary)
#     ("otheromics", "#F39B7F"), # Peach (complementary)
# ]

# # Okabe-Ito (Colorblind Safe & High Contrast)
# DEFAULT_COLOR_RULES = [
#     ("rna-seq", "#E69F00"),  # Orange (high contrast)
#     ("atac-seq", "#56B4E9"), # Sky blue
#     ("chip-seq", "#009E73"), # Blue-green
#     ("derived", "#F0E442"),  # Bright yellow (note: requires white background)
#     ("otheromics", "#CC79A7"), # Purple-red
# ]

DEFAULT_COLOR_RULES = [
    ("rna-seq", "#E64B35"),  # Brick red (emphasis)
    ("atac-seq", "#56B4E9"), # Sky blue
    ("chip-seq", "#009E73"), # Blue-green
    ("derived", "#F0E442"),  # Bright yellow (note: requires white background)
    ("otheromics", "#CC79A7"), # Purple-red
]

ANNOTATION_KEYWORDS = ("annotation",)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply category-based colors to JBrowse tracks")
    parser.add_argument(
        "--config",
        default="config.json",
        help="Path to JBrowse config.json (default: config.json)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show prospective changes without writing config.json",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-track updates",
    )
    return parser.parse_args()


def load_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"config.json not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def should_skip_annotation(category_value: str) -> bool:
    lowered = category_value.lower()
    return any(key in lowered for key in ANNOTATION_KEYWORDS)


def pick_color(category_value: str) -> Optional[str]:
    lowered = category_value.lower()
    if should_skip_annotation(lowered):
        return None
    for keyword, color in DEFAULT_COLOR_RULES:
        if keyword in lowered:
            return color
    # Fall back to "Other" palette if nothing matched
    return "#74C476"


def ensure_display(track: Dict[str, Any], display_type: str, default_renderer: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    displays = track.setdefault("displays", [])
    for display in displays:
        if display.get("type") == display_type:
            return display
    new_display: Dict[str, Any] = {"type": display_type}
    if default_renderer:
        new_display.update(default_renderer)
    displays.append(new_display)
    return new_display


def _assign_renderer_poscolor(renderer: Dict[str, Any], color: str) -> bool:
    changed = False
    if renderer.get("posColor") != color:
        renderer["posColor"] = color
        changed = True
    if renderer.get("color"):
        renderer.pop("color")
        changed = True
    if renderer.get("pos_color"):
        renderer.pop("pos_color")
        changed = True
    return changed


def apply_color_to_track(track: Dict[str, Any], color: str) -> bool:
    track_type = track.get("type")
    if track_type != "QuantitativeTrack":
        return False

    changed = False
    display = ensure_display(
        track,
        "LinearWiggleDisplay",
        {
            "renderers": {
                "XYPlotRenderer": {
                    "type": "XYPlotRenderer",
                }
            }
        },
    )

    renderers = display.setdefault(
        "renderers",
        {"XYPlotRenderer": {"type": "XYPlotRenderer"}},
    )
    for renderer in list(renderers.values()):
        if isinstance(renderer, dict):
            if _assign_renderer_poscolor(renderer, color):
                changed = True

    for field in ("color", "posColor", "pos_color"):
        if display.get(field):
            display.pop(field)
            changed = True

    return changed


def extract_category(track: Dict[str, Any]) -> Optional[str]:
    metadata = track.get("metadata") or {}
    category = metadata.get("category")
    if category:
        return category
    category_list = track.get("category")
    if isinstance(category_list, list) and category_list:
        return category_list[0]
    return None


def main() -> None:
    args = parse_args()
    config_path = Path(args.config)
    config = load_config(config_path)
    tracks = config.get("tracks", [])

    total = len(tracks)
    updated = 0
    skipped = 0

    for track in tracks:
        category_value = extract_category(track)
        if not category_value:
            skipped += 1
            continue
        color = pick_color(category_value)
        if color is None:
            skipped += 1
            continue

        if apply_color_to_track(track, color):
            updated += 1
            if args.verbose:
                print(f"[updated] {track.get('trackId')} -> {category_value} ({color})")
        else:
            skipped += 1

    print(f"Processed tracks: {total}, updated: {updated}, skipped: {skipped}")

    if args.dry_run:
        print("Dry run enabled; config.json not modified.")
        return

    tmp_fd, tmp_name = tempfile.mkstemp(prefix=f"{config_path.name}.", suffix=".tmp", dir=str(config_path.parent))
    try:
        with os.fdopen(tmp_fd, "w", encoding="utf-8") as fh:
            json.dump(config, fh, indent=2)
            fh.write("\n")
        os.replace(tmp_name, config_path)
    except Exception:
        Path(tmp_name).unlink(missing_ok=True)
        raise
    print(f"Saved updated config to {config_path}")


if __name__ == "__main__":
    main()
