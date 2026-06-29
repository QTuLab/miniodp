# Portal Demo

This directory defines the reference layout for the portal demo package.
The package is prepared from a zebrafish deployment and distributed as a separate download.

## Scope

This demo package contains:
- one species: `zebrafish`
- Gene Info
- Bulk RNA
- scRNA
- scATAC
- BulkMulti
- JBrowse 2
- SequenceServer

## Download

- file: `portal_demo.tar.gz`
- url: `https://tulab.genetics.ac.cn/~qtu/miniodp/portal_demo.tar.gz`
- size: `1.7G`
- SHA256: `85f5d2502744bc740602aeaecc3e0d9b5665c080fa5bab923d549107377d223c`

The checksum list is published at:

- `https://tulab.genetics.ac.cn/~qtu/miniodp/SHA256SUMS.txt`

## Layout

```text
demo/portal_demo/
├── dash/
│   └── data/
│       ├── species_adapters.toml
│       └── zebrafish/
│           ├── geneinfo.db
│           ├── geneinfo.toml
│           ├── bulkRNA/
│           ├── BulkMulti/
│           ├── scRNA/
│           └── scATAC/
├── hugo/
│   └── data/
│       └── species_display.toml
├── jbrowse2/
│   └── data/
│       └── Danio_rerio/
└── sequenceserver/
    └── data/
```

## Selected datasets

- scRNA:
  - `2023_DevCell_Farrell/003_004hpf`
- scATAC:
  - `2022_NatCommun_Crump/cranial_nc_005dpf`
- BulkMulti:
  - `2025Bulk`
- JBrowse 2:
  - `Danio_rerio`
- SequenceServer:
  - `Danio_rerio.GRCz11.cdna.all.fa*`

## Intended use

Use this subset to validate:
- species discovery in Dash
- gene search with `geneinfo.db`
- Bulk RNA loading from Parquet
- scRNA sample discovery and plotting
- scATAC sample discovery and plotting
- BulkMulti discovery and region-linked visualization
- JBrowse 2 loading for the zebrafish demo genome
- SequenceServer startup with a zebrafish cDNA database

## Integration notes

To test against the current repository:

1. Download and verify the package from the repository root:

```bash
curl -L -o portal_demo.tar.gz https://tulab.genetics.ac.cn/~qtu/miniodp/portal_demo.tar.gz
curl -L -o SHA256SUMS.txt https://tulab.genetics.ac.cn/~qtu/miniodp/SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt --ignore-missing
tar -xzf portal_demo.tar.gz -C demo
```

2. For a local Dash test without copying data into the repository:

```bash
cd dash
DASH_DATA_PATH=../demo/portal_demo/dash/data python app.py
```

3. To use Docker Compose, copy or sync the demo Dash data into `dash/data/`
   because the compose file mounts that directory into the container.

4. To build a single-species Hugo demo, copy the demo
   `species_display.toml` over `hugo/data/species_display.toml` in a temporary
   working tree, then build with:

```bash
cd hugo
hugo --config config/hugo_default.toml --baseURL http://localhost:1313/miniodp/
```

5. Copy `demo/portal_demo/sequenceserver/data/` into `sequenceserver/data/`
   before starting SequenceServer.

6. If you want to test the repository JBrowse helper script directly against
   the demo bundle, run:

```bash
python jbrowse2/scripts/generate_loader.py Danio_rerio \
  --data-root demo/portal_demo/jbrowse2/data \
  --app-root demo/portal_demo/jbrowse2
```

This command is for repository-local validation. On a deployed JBrowse 2 host, keep using the real app and data roots.
