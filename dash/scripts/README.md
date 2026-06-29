# Dash Scripts (Entry Points and Style Guide)

This directory contains scripts for Dash runtime and data preparation. The repository has unified on the **h5ad** data system; historical migration scripts have been moved to `dash/archive/` (no longer maintained as entry points).

## Supported Scripts

### 1) scRNA: Seurat → h5ad

- Script: `dash/scripts/scrna_convert_h5ad.py`
- Functionality: Read Seurat objects from `.Rdata/.rds` and generate one `.h5ad` per sample; supports incremental MAGIC version generation; optionally generates feather cache for UMAP+obs.
- Usage examples:
  - `python dash/scripts/scrna_convert_h5ad.py --help`
  - `python dash/scripts/scrna_convert_h5ad.py convert ...`
  - `python dash/scripts/scrna_convert_h5ad.py magic ...`

### 2) scATAC: Gene Activity → h5ad

- Script: `dash/scripts/scatac_convert_h5ad.py`
- Functionality: Extract Gene Activity (log1p) from Signac/Seurat objects, generate one `.h5ad` + `*_metadata.json` per sample; strict validation of data scale to avoid using counts/scaled matrices as input.
- Usage examples:
  - `python dash/scripts/scatac_convert_h5ad.py --help`

### 3) Single-cell Statistics: Update `species_display.toml`

- Script: `dash/scripts/summarize_single_cell.py`
- Functionality: Scan sample metadata under `dash/data/<species>/{scRNA,scATAC}` (skipping studies whose directory name contains `_Unpub_`), aggregate sample and cell counts, write to `dash/data/<species>/single_cell_stats.csv`, and update statistics fields for the corresponding species in `dash/data/species_display.toml`. By default, the script also mirrors the updated totals into the repo copy of `hugo/data/species_display.toml` for the next Hugo build.
- Usage examples:
  - `python dash/scripts/summarize_single_cell.py --help`

## Other Scripts/Subdirectories

- `dash/scripts/build_geneinfo_from_local.py`: Build `geneinfo.db` + `geneinfo.toml` from local GTF, optional gene-level annotation TSV/ZIP, and optional GO/KEGG TSV files. Useful for species that do not have an Ensembl-ready source but already have a stable local annotation set.
- `dash/scripts/bulkrna_convert_from_rsem.py`: Convert RSEM `*.genes.results` files to Dash `bulkRNA` Parquet files.
- `dash/scripts/geneinfo_from_ensembl/`: Build `geneinfo.db` + `geneinfo.toml` online from Ensembl (requires network access).
- `dash/scripts/igdb_name_overrides/`: Create and import IGDB name override table `IGDB_NAME` for Medaka (manually maintained TSV).
- `dash/scripts/landscape_view/`: Landscape View offline rendering CLI and config template.
- `dash/scripts/logging_utils.py`: Unified logging utility for scripts (print proxy).

## Style Guide (for future script additions)

### Code Language and Output

- **Code comments, logs, CLI help, and error messages should consistently use English** (aligned with repository standards) to avoid maintenance costs from mixed languages.
- Documentation (`README.md` / `docs/*.md`) should be in English for public repository; commands/paths remain unchanged.

### Paths and Configuration

- Use `pathlib.Path` for path handling; avoid hardcoding old directory structures (e.g., `services/...`).
- Data root directory follows `DASH_DATA_PATH` (default `dash/data`), consistent between scripts and services.

### Logging and Maintainability

- Script entry points should use `#!/usr/bin/env python3` + `argparse`, must support `--help`.
- Use `dash/scripts/logging_utils.py` consistently:
  - `logger = configure_script_logger(__name__)`
  - `print = create_print_proxy(logger)`
- Strict failure: Error out directly when input is missing/structure is inconsistent, do not silently fall back to generating "seemingly usable" false results.
