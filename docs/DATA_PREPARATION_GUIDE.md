# Data Preparation Guide

## Objectives
- Standardize datasets for Dash to query quickly.
- Keep species-specific details in configuration, not in code.
- Keep public statistics traceable to explicit run, base, sample-group, and
  cell counts.

For the end-to-end sequence of adding a new species, see
[ONBOARDING_GUIDE.md](ONBOARDING_GUIDE.md).

## Data directory layout
- `dash/data/{species}/geneinfo.db`
- `dash/data/{species}/geneinfo.toml`
- `dash/data/{species}/bulkRNA/`
- `dash/data/{species}/scRNA/`
- `dash/data/{species}/scATAC/`
- `dash/data/{species}/BulkMulti/`
- `dash/data/{species}/misc/` (optional auxiliary files)

## Supported formats
- SQLite for tabular metadata and gene info
- Parquet for large tables where columnar reads help
- H5AD for single-cell matrices
- BigWig for signal tracks

## Preparation workflow (outline)
1) Collect raw inputs per species (expression matrices, metadata, peaks, annotations).
2) Use `dash/scripts/` converters to generate SQLite/Parquet/H5AD/BigWig outputs.
3) Validate schema and basic stats; keep logs with command lines used.
4) Place outputs under `dash/data/{species}/` following the layout above.
5) Update display/adaptor config as needed (`dash/data/species_adapters.toml`, `hugo/data/species_display.toml`).

## Raw Input Checklist

Use this checklist before conversion:

- Genome FASTA and GTF/GFF3 annotation.
- Gene-level annotation TSV with gene ID, gene symbol, and description.
- Optional functional annotation TSV files for GO, KEGG, InterPro,
  transcription factor family, and orthologs.
- Bulk RNA quantification output, such as StringTie `*_gene_abundance.tab` or
  RSEM `*.genes.results`.
- Bulk RNA metadata with `sample_id`, `source`, `sample_name`, and optional
  `description`.
- ATAC-seq and ChIP-seq BigWig and peak files for JBrowse.
- Track metadata in `MetaData.csv` with stable labels and categories.
- Single-cell Seurat, Signac, or H5AD data with cell metadata, gene IDs, UMAP,
  and cell-type labels.
- BLAST FASTA files for genome, cDNA, and protein searches.

## Core Dash Outputs

### Gene Info

Build from local annotation:

```bash
python dash/scripts/build_geneinfo_from_local.py \
  --gtf annotation.gtf \
  --annotation-tsv gene_annotation.tsv \
  --download-go-obo \
  --species-label "Genus species" \
  -o dash/data/{species}
```

Expected outputs:

- `geneinfo.db`
- `geneinfo.toml`

### Bulk RNA

For StringTie:

```bash
python dash/scripts/bulkrna_convert_from_stringtie.py \
  --input-dir quant \
  --metadata bulkrna_metadata.csv \
  --output-dir dash/data/{species}/bulkRNA \
  --species {species} \
  --validate-genes
```

For RSEM, use `dash/scripts/bulkrna_convert_from_rsem.py`.

### scRNA

```bash
python dash/scripts/scrna_convert_h5ad.py convert \
  --input-dir seurat_objects \
  --output-dir dash/data/{species}/scRNA \
  --study 2024_Journal_Author
```

### scATAC

```bash
python dash/scripts/scatac_convert_h5ad.py \
  --input signac_object.rds \
  --output-dir dash/data/{species}/scATAC \
  --study 2024_Journal_Author \
  --assay ACTIVITY \
  --slot data
```

### BulkMulti

```bash
python dash/scripts/bulkmulti_convert_from_pipeline.py \
  --species {species} \
  --pipeline-dir pipeline \
  --metadata bulkmulti_metadata.csv \
  --gtf annotation.gtf \
  --output-dir dash/data/{species}/BulkMulti
```

## Public Statistics

The Hugo home page uses four public statistics:

- `bulk_datasets`: number of bulk sequencing runs.
- `bulk_bases`: sequenced bases from bulk assays.
- `single_cell_datasets`: number of single-cell sample groups.
- `cells`: number of cells in published single-cell datasets.

Keep run-level audit tables outside the public repository if they contain
unpublished samples. When values are updated, record the exact source table,
command, and date in your local maintenance notes.

## Quality checks
- Spot-check row counts and key columns in SQLite/Parquet.
- Ensure gene identifiers match between modalities.
- Confirm BigWig chromosome names align with browser expectations.
- Confirm `species_display.toml` links open the expected Dash, JBrowse, and
  SequenceServer routes.
- Run converter commands once with `--dry-run` when the script supports it.
