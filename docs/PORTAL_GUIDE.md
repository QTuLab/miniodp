# Portal Guide

This guide explains how the public-facing `miniodp` portal is assembled from the repository components.

## Components

The portal stack has four user-facing parts and one shared proxy layer:

- `hugo/`: static home page and species pages
- `dash/`: interactive data applications
- `jbrowse2/`: genome browser app and track resources
- `sequenceserver/`: BLAST search service
- `common/nginx/miniodp.conf`: reverse-proxy and static route reference

## Recommended deployment model

The current deployment model is hybrid:

- `Nginx` runs on the host
- `Hugo` builds static files on the host
- `JBrowse 2` is served as static content on the host
- `Dash` runs in Docker
- `SequenceServer` runs in Docker

This layout keeps the static parts simple and isolates the long-running application services.

## Component relationships

### Hugo

`hugo/` provides the public entry pages and species landing pages.

Typical build command:

```bash
cd hugo
hugo --minify --config config/hugo_default.toml --baseURL https://example.org/miniodp/
```

For local preview, use:

```bash
cd hugo
hugo server --config config/hugo_default.toml --baseURL http://localhost:1313/miniodp/ --port 1313
```

The generated site is expected to be served behind `/miniodp/`.
Private deployment override files are not included in this public repository.
Use `--baseURL` or a local untracked config file for environment-specific
values.

### Dash

`dash/` provides the interactive applications for:

- gene information
- bulk RNA
- scRNA
- scATAC
- BulkMulti

Typical startup:

```bash
cd dash
docker compose up -d
```

The compose file binds the service to `127.0.0.1:3939` and expects runtime data under `dash/data/`.

### JBrowse 2

`jbrowse2/` contains helper scripts used to prepare reference assemblies and tracks for a separate JBrowse 2 app directory.

Typical responsibilities:

- prepare the reference sequence layout
- generate loader scripts from metadata
- register tracks in the app
- build text indexes

The final static app is expected under `/miniodp/jbrowse2/`.

### SequenceServer

`sequenceserver/` provides the BLAST service entry point.

Typical startup:

```bash
cd sequenceserver
docker compose up -d
```

The compose file uses the official `wurmlab/sequenceserver` image, binds to `127.0.0.1:4567`, and mounts `sequenceserver/data/` into `/db`.

## Nginx routing

The repository includes a reference configuration snippet:

```text
common/nginx/miniodp.conf
```

It defines routes for:

- `/miniodp/` to Hugo static files
- `/miniodp/jbrowse2/` to JBrowse 2 static resources
- `/miniodp/dash/` to the Dash container
- `/miniodp/sequenceserver/` to the SequenceServer container

## Portal demo package

The portal demo is distributed separately from the repository and is intended to validate the wiring between all portal components.

Current download:

- `https://tulab.genetics.ac.cn/~qtu/miniodp/portal_demo.tar.gz`
- SHA256: `85f5d2502744bc740602aeaecc3e0d9b5665c080fa5bab923d549107377d223c`

Its documented layout is:

```text
portal_demo/
├── dash/
│   └── data/
├── hugo/
│   └── data/
├── jbrowse2/
│   └── data/
└── sequenceserver/
    └── data/
```

See:

- [../demo/portal_demo/README.md](../demo/portal_demo/README.md)

## Demo validation path

To validate the portal with the demo package:

1. Build Hugo with the demo `species_display.toml`.
2. Point Dash to the demo `dash/data/` directory.
3. Copy the demo JBrowse 2 species bundle into the JBrowse 2 app data directory and register the assembly and tracks.
4. Copy the demo BLAST database files into `sequenceserver/data/`.
5. Start Dash and SequenceServer with Docker Compose.
6. Serve the static routes and reverse proxies with Nginx.

## Related documents

- [../README.md](../README.md)
- [../dash/README.md](../dash/README.md)
- [../hugo/README.md](../hugo/README.md)
- [../jbrowse2/README.md](../jbrowse2/README.md)
- [ONBOARDING_GUIDE.md](ONBOARDING_GUIDE.md)
- [DATA_PREPARATION_GUIDE.md](DATA_PREPARATION_GUIDE.md)
- [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md)
