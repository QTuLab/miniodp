# SequenceServer

This directory provides the Docker entry point for the BLAST service used by miniodp.

## Runtime model

- Container image: official `wurmlab/sequenceserver`
- Mounted database path inside the container: `/db`
- Public path behind Nginx: `/miniodp/sequenceserver/`

## Local startup

Populate `sequenceserver/data/` with BLAST databases, then start:

```bash
python3 scripts/manage_blastdb.py materialize
bash scripts/build_all_blastdbs.sh
docker compose up -d
```

`scripts/manage_blastdb.py` reads `config/blastdb_manifest.csv` and can:

- download managed FASTA files from official sources
- copy locally prepared FASTA files when a manifest row uses `source_kind=local`
- build BLAST indexes for files present in `sequenceserver/data/`

Useful examples:

```bash
# Preview managed records
python3 scripts/manage_blastdb.py list

# Refresh only managed FASTA files for selected species
python3 scripts/manage_blastdb.py materialize --species zebrafish --species medaka --species mexican_tetra --species cattle --species hydra --species lancelet

# Rebuild databases for all FASTA files currently present
bash scripts/build_all_blastdbs.sh
```

The build step skips FASTA files that are not present in `sequenceserver/data/`.
Unmanaged entries, such as locally prepared genome bundles, can be copied into
`sequenceserver/data/` manually before running the build step.

Runtime configuration files are kept in `config/`:

- `blastdb_manifest.csv`: source FASTA files and BLAST database build metadata
- `jbrowse_genome_links.csv`: genome BLAST database to JBrowse assembly mapping

## JBrowse links for genome hits

SequenceServer loads `scripts/jbrowse_links.rb` at startup through the
container path `sequenceserver -r /opt/sequenceserver/jbrowse_links.rb`.
The extension adds a `JBrowse` link to BLAST hits that come from genome
databases. It does not add links for cDNA or protein hits, because their hit
coordinates are transcript or protein coordinates rather than genomic
coordinates.

The mapping from BLAST genome databases to JBrowse assemblies is maintained in:

```text
config/jbrowse_genome_links.csv
```

Each enabled row records the BLAST FASTA file name, the database title, the
target JBrowse assembly, default annotation tracks, and optional reference-name
prefixes to strip when BLAST sequence IDs differ from JBrowse reference names.
During result rendering, the extension asks SequenceServer which BLAST database
contains the hit, looks up the matching row, and builds a root-relative URL such
as:

```text
/miniodp/jbrowse2/?assembly=Danio_rerio&loc=14:26709229..26719269
```

The reference name in `loc` should match the JBrowse reference FASTA. The
extension handles common FASTA header patterns from Ensembl and NCBI, including
`dna_rm:chromosome`, `dna_rm:primary_assembly`, and `refseq_seqid=NC_...`.
If a new genome uses another header convention, add a regression case to
`scripts/check_jbrowse_links.rb` before deployment.

The generated JBrowse view opens the configured annotation tracks plus a dynamic
`BLAST HSPs` session track. The default location window is about twice the
combined HSP span, with at least `MINIODP_BLAST_JBROWSE_PADDING` bp on each side
for short hits.

The public path prefix is controlled by `MINIODP_BASE_PATH`, which defaults to
`/miniodp`. This keeps links independent of the deployment domain. If the portal
is deployed under a different path prefix, set `MINIODP_BASE_PATH` in the
SequenceServer environment and restart the container.

Lightweight validation without real BLAST databases:

```bash
ruby scripts/check_jbrowse_links.rb
```

Deployment only requires syncing the SequenceServer config files and restarting
the SequenceServer container. BLAST databases do not need to be rebuilt for this
link extension.

Files to sync for this feature:

```bash
rsync -avhP config/jbrowse_genome_links.csv <server>:/opt/miniodp/sequenceserver/config/
rsync -avhP scripts/jbrowse_links.rb scripts/check_jbrowse_links.rb <server>:/opt/miniodp/sequenceserver/scripts/
```

After syncing, restart only the SequenceServer service:

```bash
cd /opt/miniodp/sequenceserver
docker compose restart sequenceserver
```

This assumes the server already uses a compose file that mounts
`scripts/jbrowse_links.rb` and `config/jbrowse_genome_links.csv`. If that
compose structure has not been deployed yet, update `docker-compose.yml` first
and start the service once with the normal deployment command. Later changes to
the Ruby or CSV files only require `docker compose restart sequenceserver`.

## Portal demo

The public portal demo package ships a zebrafish cDNA BLAST database under:

```text
portal_demo/
└── sequenceserver/
    └── data/
        └── Danio_rerio.GRCz11.cdna.all.fa*
```

Copy those files into `sequenceserver/data/` before starting the service.
