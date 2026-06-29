# JBrowse 2 Preparation

`jbrowse2/` contains helper scripts for preparing reference assemblies and
tracks for a JBrowse 2 static app. The app itself is expected to be installed or
served outside this repository.

## Expected Species Bundle

Each species should be prepared under the JBrowse app data root:

```text
data/<assembly_name>/
├── MetaData.csv
├── reference/
│   └── genome.fa
└── tracks/
    ├── annotation/
    ├── rnaseq/
    ├── atacseq/
    ├── chipseq/
    ├── derived/
    └── otheromics/
```

`MetaData.csv` is the authoritative track list used by
`scripts/generate_loader.py`. Keep labels stable because Dash and Hugo links may
refer to those track names.

## Prepare Reference

```bash
python jbrowse2/scripts/prepare_reference.py <assembly_name> \
  --data-root /path/to/jbrowse2/data \
  --app-root /path/to/jbrowse2 \
  --assembly-name <assembly_name>
```

This creates FASTA indexes and registers the assembly in the JBrowse
`config.json`.

## Prepare Tracks

```bash
python jbrowse2/scripts/prepare_tracks.py <assembly_name> \
  --data-root /path/to/jbrowse2/data
```

This converts GFF3/GFF files to tabix-indexed files and BED files to BigBed
when the required command-line tools are available.

## Register Tracks

```bash
python jbrowse2/scripts/generate_loader.py <assembly_name> \
  --data-root /path/to/jbrowse2/data \
  --app-root /path/to/jbrowse2

bash /path/to/jbrowse2/data/<assembly_name>/run_add_tracks.sh
```

The generated script skips tracks that are already present in `config.json`.

## Maintenance Helpers

- `scripts/count_tracks.py`: count tracks by assembly and category.
- `scripts/reorder_tracks_by_metadata.py`: reorder existing tracks to match
  `MetaData.csv`.
- `scripts/auto_color_tracks.py`: apply category-based colors to quantitative
  tracks.

When production `config.json` is maintained on a server, copy it locally first,
run these scripts against the local copy, review the diff, and then sync the
updated `config.json` with the new track files.

