# Onboarding Guide

This guide is written for a new group that wants to understand, test, and
extend `miniodp` without prior knowledge of the hosted QTu Lab deployment.

## What You Need to Provide

A new species module usually needs these inputs before it can become a useful
portal page:

- Basic species profile: display name, scientific name, short description, one
  representative image, and the public URL prefix used by your deployment.
- Genome reference: FASTA, chromosome sizes, and GFF3 or GTF annotation.
- Gene annotation table: gene identifier, gene symbol, description, and
  optional functional annotations such as GO, KEGG, InterPro, transcription
  factor family, and ortholog links.
- Bulk RNA-seq data: one quantification table per sample group, metadata with
  study/source, sample name, description, accessions, run list, and read/base
  statistics.
- ATAC-seq and ChIP-seq data: BigWig signal files, peak files, and metadata for
  study/source, sample group, assay, antibody or mark when relevant, accessions,
  run list, and read/base statistics.
- Single-cell data: Seurat object, H5AD file, or equivalent matrix with cell
  metadata, UMAP coordinates, gene identifiers, sample group labels, cell-type
  labels, and study metadata.
- BLAST data: genome, cDNA, and protein FASTA files when available.

Large runtime datasets are not stored in this repository. The repository
contains the code, configuration templates, and conversion scripts needed to
build those runtime datasets.

## Recommended First Test

Before adding your own species, test the public demo package.

```bash
curl -L -o portal_demo.tar.gz https://tulab.genetics.ac.cn/~qtu/miniodp/portal_demo.tar.gz
curl -L -o SHA256SUMS.txt https://tulab.genetics.ac.cn/~qtu/miniodp/SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt --ignore-missing
tar -xzf portal_demo.tar.gz -C demo
```

Then follow the component-specific notes in:

- [../demo/portal_demo/README.md](../demo/portal_demo/README.md)
- [PORTAL_GUIDE.md](PORTAL_GUIDE.md)
- [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md)

The demo package is a validation bundle. It is intended to prove that the
components can find their data, not to replace a full production deployment.

## Naming Decisions

Choose stable names before preparing data.

- `species_key`: lower-case Dash and Hugo key, such as `zebrafish`.
- `assembly_name`: JBrowse assembly name, such as `Danio_rerio`.
- `study_id`: stable public identifier for a dataset group. Prefer
  `Year_Journal_CorrespondingAuthor` and add a fourth segment only when one
  paper has multiple independent data groups, for example
  `2020_Science_Brunet_ATAC`.
- `sample_id`: machine-readable sample group identifier used by converted
  files.
- `sample_name`: short human-readable label shown in Dash.

For studies with more than one BioProject or clearly independent assay block,
keep each BioProject or assay block as a separate `study_id` when that makes
data processing and provenance easier to audit.

## Add a Species Page

Edit [../hugo/data/species_display.toml](../hugo/data/species_display.toml).
Add one table for your species:

```toml
[new_species]
name = "Common name"
scientific_name = "Genus species"
image = "new_species.jpg"
description = "Short public description."
jbrowse_url = "/miniodp/jbrowse2/?assembly=Genus_species"
dash_url = "/miniodp/dash/new_species"
blast_url = "/miniodp/sequenceserver/"
order = 99
bulk_datasets = 0
bulk_bases = 0
single_cell_datasets = 0
cells = 0
```

Place the image under `hugo/static/images/`. Build locally with:

```bash
cd hugo
hugo --config config/hugo_default.toml --baseURL http://localhost:1313/miniodp/
```

For deployment, set `--baseURL` to your public portal URL.

## Add Dash Gene Search

Build the Dash gene information database from a local annotation set:

```bash
python dash/scripts/build_geneinfo_from_local.py \
  --gtf /path/to/annotation.gtf \
  --annotation-tsv /path/to/gene_annotation.tsv \
  --go-kegg-tsv /path/to/gene2go_kegg.tsv \
  --download-go-obo \
  --species-label "Genus species" \
  -o dash/data/new_species
```

The output directory must contain:

```text
dash/data/new_species/
├── geneinfo.db
└── geneinfo.toml
```

If your gene identifiers are not standard Ensembl-like identifiers, configure
the adapter in [../dash/data/species_adapters.toml](../dash/data/species_adapters.toml).
Use `gene_id` for generic numeric or locus identifiers, and use `dual_system`
only when one species intentionally has two linked gene ID systems.

## Add Bulk RNA Data

Prepare one metadata CSV using
[../dash/scripts/bulkrna_metadata_template.csv](../dash/scripts/bulkrna_metadata_template.csv).
Required columns are:

- `sample_id`: must match the quantification file prefix.
- `source`: study or dataset label shown in Dash.
- `sample_name`: sample group label shown in Dash.
- `description`: optional sample description.

For StringTie output:

```bash
python dash/scripts/bulkrna_convert_from_stringtie.py \
  --input-dir /path/to/stringtie_quant \
  --metadata /path/to/bulkrna_metadata.csv \
  --output-dir dash/data/new_species/bulkRNA \
  --species new_species \
  --validate-genes
```

For RSEM output, use `dash/scripts/bulkrna_convert_from_rsem.py` with the same
metadata idea.

## Add Single-Cell Data

For scRNA Seurat objects:

```bash
python dash/scripts/scrna_convert_h5ad.py convert \
  --input-dir /path/to/seurat_objects \
  --output-dir dash/data/new_species/scRNA \
  --study 2024_Journal_Author
```

For scATAC gene-activity data:

```bash
python dash/scripts/scatac_convert_h5ad.py \
  --input /path/to/signac_or_seurat_object.rds \
  --output-dir dash/data/new_species/scATAC \
  --study 2024_Journal_Author \
  --assay ACTIVITY \
  --slot data
```

After adding single-cell data, update the public statistics:

```bash
cp hugo/data/species_display.toml dash/data/species_display.toml
python dash/scripts/summarize_single_cell.py --species new_species
```

The script updates `dash/data/species_display.toml` and mirrors the totals back
to `hugo/data/species_display.toml` by default.

## Add BulkMulti Data

BulkMulti links processed RNA, ATAC, ChIP, enhancer-like signature, and gene
regulatory network outputs in the landscape view.

Prepare metadata from
[../dash/scripts/bulkmulti_metadata_template.csv](../dash/scripts/bulkmulti_metadata_template.csv),
then run:

```bash
python dash/scripts/bulkmulti_convert_from_pipeline.py \
  --species new_species \
  --pipeline-dir /path/to/pipeline \
  --metadata /path/to/bulkmulti_metadata.csv \
  --gtf /path/to/annotation.gtf \
  --output-dir dash/data/new_species/BulkMulti
```

Use explicit path override columns in the metadata file when your pipeline
outputs do not follow the default naming convention.

## Add JBrowse Reference and Tracks

Prepare a species bundle under a JBrowse app data root:

```text
data/Genus_species/
├── MetaData.csv
├── reference/
│   └── genome.fa
└── tracks/
    ├── annotation/
    ├── rnaseq/
    ├── atacseq/
    ├── chipseq/
    └── derived/
```

Run the helper scripts from the repository root:

```bash
python jbrowse2/scripts/prepare_reference.py Genus_species \
  --data-root /path/to/jbrowse2/data \
  --app-root /path/to/jbrowse2 \
  --assembly-name Genus_species

python jbrowse2/scripts/prepare_tracks.py Genus_species \
  --data-root /path/to/jbrowse2/data

python jbrowse2/scripts/generate_loader.py Genus_species \
  --data-root /path/to/jbrowse2/data \
  --app-root /path/to/jbrowse2

bash /path/to/jbrowse2/data/Genus_species/run_add_tracks.sh
```

If the generated `config.json` is managed on another host, run the scripts
against a local copy of that app and then sync the updated `config.json` and
track files together.

## Add SequenceServer Databases

Edit [../sequenceserver/config/blastdb_manifest.csv](../sequenceserver/config/blastdb_manifest.csv)
for genome, cDNA, and protein FASTA files. Managed rows can download from a
URL; local rows expect you to copy the FASTA file into `sequenceserver/data/`.

Then run:

```bash
cd sequenceserver
python3 scripts/manage_blastdb.py materialize --species new_species
bash scripts/build_all_blastdbs.sh
docker compose restart sequenceserver
```

For genome BLAST hits to open in JBrowse, also add a row to
[../sequenceserver/config/jbrowse_genome_links.csv](../sequenceserver/config/jbrowse_genome_links.csv).
Validate link generation with:

```bash
ruby sequenceserver/scripts/check_jbrowse_links.rb
```

## Final Validation Checklist

Before publishing a species module, verify:

- Hugo builds without missing-file warnings caused by custom config files.
- The species card links to the expected Dash, JBrowse, and BLAST routes.
- Dash discovers the species under `DASH_DATA_PATH`.
- Gene search returns gene name, gene ID, locus, and description.
- Bulk RNA, scRNA, scATAC, and BulkMulti pages list expected studies and samples.
- JBrowse opens the selected assembly and track names from `species_display.toml`.
- SequenceServer can load the FASTA databases and genome BLAST links open JBrowse.
- Public statistics use `bulk_datasets`, `bulk_bases`, `single_cell_datasets`,
  and `cells`; internal run/base audit tables should be kept outside this
  public repository when they contain unpublished data.
