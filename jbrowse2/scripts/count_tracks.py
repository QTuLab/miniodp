#!/usr/bin/env python3
"""
Count tracks in JBrowse 2 config.json.

Features:
- Count tracks by species (assembly)
- Support category-based statistics (--by-category)

Usage:
    python3 count_tracks.py
    python3 count_tracks.py --by-category
    python3 count_tracks.py /path/to/config.json --by-category
"""
import json
import sys
import argparse
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Count tracks per species in JBrowse 2 config.json"
    )
    parser.add_argument(
        "config_path",
        nargs="?",
        default="config.json",
        help="Path to config.json (default: config.json)",
    )
    parser.add_argument(
        "--by-category",
        "-c",
        action="store_true",
        help="Statistics by category",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    path = Path(args.config_path)

    if not path.exists():
        print(f"Error: File not found -> {path}")
        sys.exit(1)

    try:
        print(f"Reading config: {path} ...")
        with open(path, "r", encoding="utf-8") as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error: JSON parse failed -> {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    registered_assemblies = set()
    if "assemblies" in config:
        for asm in config["assemblies"]:
            if "name" in asm:
                registered_assemblies.add(asm["name"])

    print(f"Registered assemblies: {len(registered_assemblies)}")

    tracks = config.get("tracks", [])
    total_tracks = len(tracks)

    if args.by_category:
        print_by_category(tracks, registered_assemblies)
    else:
        print_by_assembly(tracks, registered_assemblies)

    print(f"\nTotal tracks: {total_tracks}")


def print_by_assembly(tracks, registered_assemblies):
    """Statistics by species"""
    stats = defaultdict(int)
    for name in registered_assemblies:
        stats[name] = 0

    for track in tracks:
        assembly_names = track.get("assemblyNames", [])
        for asm_name in assembly_names:
            stats[asm_name] += 1

    print("-" * 50)
    print(f"{'Assembly Name':<30} | {'Track Count':<10}")
    print("-" * 50)

    sorted_stats = sorted(stats.items(), key=lambda item: (-item[1], item[0]))

    for asm_name, count in sorted_stats:
        flag = ""
        if asm_name not in registered_assemblies:
            flag = "(not defined in assemblies)"
        print(f"{asm_name:<30} | {count:<10} {flag}")

    print("-" * 50)
    print(f"Total references: {sum(stats.values())}")


def print_by_category(tracks, registered_assemblies):
    """Statistics by species and category"""
    # stats[assembly][category] = count
    stats = defaultdict(lambda: defaultdict(int))

    for track in tracks:
        assembly_names = track.get("assemblyNames", [])
        category = track.get("category", ["Uncategorized"])
        if isinstance(category, list):
            category = category[0] if category else "Uncategorized"
        for asm_name in assembly_names:
            stats[asm_name][category] += 1

    sorted_assemblies = sorted(stats.keys())

    for asm_name in sorted_assemblies:
        flag = ""
        if asm_name not in registered_assemblies:
            flag = " (not defined in assemblies)"

        categories = stats[asm_name]
        total = sum(categories.values())

        print("")
        print("=" * 50)
        print(f"Assembly: {asm_name}{flag}")
        print("-" * 50)

        sorted_categories = sorted(categories.items(), key=lambda x: x[0])
        max_cat_len = max(len(cat) for cat, _ in sorted_categories) if sorted_categories else 15

        for cat, count in sorted_categories:
            print(f"  {cat:<{max_cat_len}} : {count:>5}")

        print("-" * 50)
        print(f"  {'Total':<{max_cat_len}} : {total:>5}")


if __name__ == "__main__":
    main()
