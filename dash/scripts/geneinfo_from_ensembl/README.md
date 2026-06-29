# Geneinfo from Ensembl

## Introduction

`build_geneinfo_from_ensembl.py` is an automated tool for fetching gene annotation information from the Ensembl database and building the `geneinfo.db` database file required by the miniodp project.

## Features

- ✅ **Fully automated**: Build database with a single command
- ✅ **Authoritative data source**: Fetch latest data directly from official Ensembl API
- ✅ **Data completeness**: Includes gene info, GO annotations, orthologs, protein domains, etc.
- ✅ **Error handling**: Includes retry mechanism and detailed error reporting
- ✅ **Progress tracking**: Real-time display of data fetch progress
- ✅ **Good compatibility**: Generated database fully compatible with existing system

## Supported Data Types

### Core Data Tables
1. **ENS_INFO** - Basic gene information (gene name, type, description)
2. **ENS_Locus** - Genomic position coordinates
3. **ENS_GO** - GO functional annotations
4. **ENS_INTERPRO** - InterPro protein domains
5. **ENS_Seq** - Gene sequence information
6. **ENS_OrthoHuman** - Human orthologs
7. **ENS_OrthoMouse** - Mouse orthologs
8. **ENS_OrthoFly** - Fly orthologs

### Auto-generated
- **geneinfo.toml** - Search configuration file
- **Database indexes** - Query performance optimization

## Install Dependencies

```bash
# Use project development environment (recommended), or ensure following dependencies:
pip install requests pandas tqdm
```

## Usage

### Basic Usage

```bash
# Build zebrafish database
python build_geneinfo_from_ensembl.py danio_rerio

# Build mouse database
python build_geneinfo_from_ensembl.py mus_musculus

# Build human database
python build_geneinfo_from_ensembl.py homo_sapiens
```

### Advanced Options

```bash
# Specify output directory
python build_geneinfo_from_ensembl.py danio_rerio -o /path/to/output

# Enable verbose logging
python build_geneinfo_from_ensembl.py danio_rerio -v

# Show help
python build_geneinfo_from_ensembl.py -h
```

## Common Species Names

| Species | Ensembl Name |
|---------|--------------|
| Zebrafish | danio_rerio |
| Medaka | oryzias_latipes |
| Human | homo_sapiens |
| Mouse | mus_musculus |
| Fruit fly | drosophila_melanogaster |
| C. elegans | caenorhabditis_elegans |
| Cattle | bos_taurus |
| Arabidopsis | arabidopsis_thaliana |

## Output Files

### geneinfo.db
SQLite database file containing all gene annotation information:

```sql
-- View gene count
SELECT COUNT(*) FROM ENS_INFO;

-- View GO annotation count
SELECT COUNT(*) FROM ENS_GO;

-- Example query
SELECT * FROM ENS_INFO WHERE GeneName LIKE '%tp53%' LIMIT 5;
```

### geneinfo.toml
Search field configuration file defining searchable field types:

```toml
[[search_fields]]
label = "Gene Name"
value = "gene_name"
db_table = "ENS_INFO"
db_column = "GeneName"
match_type = "like"
```

## Data Quality Check

```bash
# Check generated database
sqlite3 geneinfo.db "
SELECT
    'ENS_INFO' as table_name, COUNT(*) as row_count FROM ENS_INFO
UNION ALL
SELECT 'ENS_GO', COUNT(*) FROM ENS_GO
UNION ALL
SELECT 'ENS_OrthoHuman', COUNT(*) FROM ENS_OrthoHuman;
"
```

## Troubleshooting

### Common Issues

1. **Invalid species name**
   ```
   ERROR - Species xxx not found in Ensembl
   ```
   - Solution: Check species name spelling, use correct Ensembl format

2. **Network connection issues**
   ```
   WARNING - Attempt 1 failed: HTTP 503
   ```
   - Solution: Check network connection, retry later (program has auto-retry mechanism)

3. **Insufficient disk space**
   ```
   ERROR - No space left on device
   ```
   - Solution: Free up disk space, database typically requires 50-500MB

4. **Permission issues**
   ```
   ERROR - Permission denied
   ```
   - Solution: Check write permissions for output directory

### Debugging Methods

```bash
# Enable verbose logging
python build_geneinfo_from_ensembl.py danio_rerio -v

# Check Ensembl service status
curl -s "https://rest.ensembl.org/info/ping"
```

## Performance Optimization

- **Batch queries**: Program uses BioMart batch API to reduce request count
- **Retry mechanism**: Auto-retry on network failure with exponential backoff
- **Index optimization**: Auto-create database indexes to improve query speed
- **Memory management**: Process large datasets in batches to avoid memory overflow

## Integration into miniodp

After build completes, place files in correct location:

```bash
# Move to species data directory
mv geneinfo.db /path/to/miniodp/dash/data/{species}/
mv geneinfo.toml /path/to/miniodp/dash/data/{species}/

# Verify integration
cd /path/to/miniodp/dash
python app.py  # Start test
```

## Update Data

```bash
# Periodic data updates (recommended quarterly)
python build_geneinfo_from_ensembl.py danio_rerio --output-dir ./backup_$(date +%Y%m%d)

# Compare old and new data
sqlite3 old_geneinfo.db "SELECT COUNT(*) FROM ENS_INFO"
sqlite3 new_geneinfo.db "SELECT COUNT(*) FROM ENS_INFO"
```

## Contact Support

If you encounter problems:
1. Check verbose logs (`-v` option)
2. Check network connection and Ensembl service status
3. Submit issue to project repository

---

*Last updated: 2025-07-15*
