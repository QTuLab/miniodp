# Landscape View CLI

This directory contains command-line tools and configuration templates for offline generation of "Landscape View (Peak-to-Gene multi-track landscape plots)".

## Files

- `landscape_view_cli.py`: CLI entry point. Reads TOML configuration, calls `BulkMultiAPI` + `LandscapeViewAPI` to generate Plotly figures and export to `html/png/pdf`.
- `landscape_view_config.toml`: Configuration template/example (species, target genes, thresholds, output directory and formats, etc.).

## Common Usage

```bash
cd dash/scripts
python landscape_view_cli.py

# Specify config file
python landscape_view_cli.py --config /path/to/landscape_view_config.toml

# Debug mode (more verbose logging)
python landscape_view_cli.py --debug
```

## Dependencies and Prerequisites

- Dependencies: `plotly` (PNG/PDF export requires `kaleido`), `toml`
- Prerequisites: Corresponding species data in `dash/data/<species>/BulkMulti/...` (DB/BigWig etc.) must be available and accessible via `DASH_DATA_PATH` or default `dash/data`.
