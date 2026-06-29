# IGDB Name Overrides (Medaka)

This directory maintains the `IGDB_NAME` table in Medaka `geneinfo.db`: provides **manually curated/supplemented gene names** (overrides) for certain IGDB genes.

## Problem Statement (Why Needed)

Medaka's gene primary key system is more complex than "standard Ensembl species": the system uses `primary_id` internally, while Medaka may have `IGDB::ENS` composite keys. Ensembl's `GeneName` is not ideal or inconsistent in some cases, thus requiring a manually maintained name table keyed by IGDB:

- Allow gene search using "IGDB curated names" (supplement/correct ENS GeneName)
- Display combined ENS name and IGDB curated name to avoid user confusion

This is not a regular data conversion pipeline, but **a small-scale, manually maintainable override layer**; the TSV file should be treated as the single source of truth.

## File Structure

```
dash/
├── scripts/
│   └── igdb_name_overrides/
│       ├── igdb_name_create.py      # Create IGDB_NAME table (one-time)
│       ├── igdb_name_import.py      # Full import from TSV (repeatable)
│       └── README.md
└── data/medaka/
    ├── geneinfo.db              # SQLite database
    └── misc/
        └── igdb_names.tsv       # Gene name file
```

## Usage Workflow

### 1) Create table (run once only)
```bash
cd dash
python scripts/igdb_name_overrides/igdb_name_create.py
```

### 2) Prepare TSV data file
Edit `data/medaka/misc/igdb_names.tsv` file with the following format:
```tsv
IGDB	GeneName	Description
IGDB07699	amh	Anti-Müllerian hormone
IGDB08949	foxl2l	Forkhead box L2-like
IGDB10884	prm	Protamine, Ref: Tamura et al. (1994)
```

### 3) Import (run on each update; overwrite mode)
```bash
cd dash
python scripts/igdb_name_overrides/igdb_name_import.py
```

## Table Structure

**IGDB_NAME table**:
- **IGDB**: Gene ID (primary key, required)
- **GeneName**: Gene name (required)
- **Description**: Description (optional)

## TSV File Format

- **IGDB**: Gene ID (required)
- **GeneName**: Gene name (required)
- **Description**: Description (optional)

## Behavioral Conventions (Important)

1. `igdb_name_import.py` uses **overwrite mode**: each import clears `IGDB_NAME` first then writes TSV content; TSV is the single source of truth.
2. Basic validation (missing columns, duplicate IGDB, empty GeneName, etc.) is performed before import; validation failure aborts to prevent dirty data.
3. To rollback: restore TSV then re-import; if concerned about mistakes, backup `geneinfo.db` first.

## Integration with Application

This table is used by `MedakaAdapter`:

- **Search**: When searching by gene name, in addition to `ENS_INFO.GeneName`, also searches `IGDB_NAME.GeneName`, supporting "IGDB curated name" queries.
- **Display**: When IGDB curated name exists, displays as `"{ens_name} [IGDB: {igdb_name}]"`; otherwise displays only ENS name.

## Update History

- **2025-07-15**: Refactored to `IGDB_NAME` table, simplified structure, unified naming convention
