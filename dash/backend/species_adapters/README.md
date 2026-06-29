# Species Adapters

Species adapters keep gene identifier differences out of the Dash user
interface code. The adapter factory reads
`dash/data/species_adapters.toml` and selects one adapter per species.

## Adapter Types

- `standard`: Ensembl-like gene identifiers stored in the `ENS` column.
- `gene_id`: generic gene identifiers such as numeric NCBI Gene IDs or local
  locus IDs. The user interface labels these as "Gene ID" rather than
  "Ensembl ID".
- `dual_system`: linked dual identifier systems. The current implementation is
  used for medaka, where IGDB identifiers are primary for expression data and
  Ensembl identifiers remain searchable.

## Configuration Example

```toml
[species.new_species]
adapter_type = "gene_id"
id_formats = ["^[0-9]+$"]
primary_key = "ENS"
```

Only add species-specific configuration when the defaults do not match the
species. Standard Ensembl-like species usually need only an `id_formats`
override.

## Medaka Dual-ID Notes

Medaka uses a dual IGDB and Ensembl model. Gene names are displayed as
`ens name [IGDB: igdb name]` when both names are available. The optional
`IGDB_NAME` table is created by:

```bash
python dash/scripts/igdb_name_overrides/igdb_name_create.py
python dash/scripts/igdb_name_overrides/igdb_name_import.py
```

This behavior is intentionally limited to the medaka adapter.

