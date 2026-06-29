# Deployment Guide

## Goals
- Provide a reproducible way to deploy the public miniodp portal.
- Keep static components simple on the host, while isolating long-running app services with containers.

## Prerequisites
- Hugo for building the static site
- Docker or Docker Compose for Dash and SequenceServer
- Nginx for reverse proxy and static hosting

## Recommended deployment model

- `Nginx` runs on the host.
- `Hugo` builds static files; it is not a long-running container service.
- `JBrowse 2` is served as static content by Nginx.
- `Dash` runs in Docker.
- `SequenceServer` runs in Docker.

This repository therefore uses a hybrid model:
- host-managed static content for `Hugo` and `JBrowse 2`
- containerized runtime services for `Dash` and `SequenceServer`

## Reference layout

- Hugo static output: `/opt/miniodp/hugo/public/`
- JBrowse 2 static resources: `/opt/miniodp/jbrowse2/`
- Dash container: `127.0.0.1:3939`
- SequenceServer container: `127.0.0.1:4567`

## Build static content

### Hugo

Build the site with the default Hugo config:

```bash
cd hugo
hugo --config config/hugo_default.toml
```

Copy the generated files to the directory served by Nginx, for example:

```bash
rsync -av --delete public/ /opt/miniodp/hugo/public/
```

### JBrowse 2

Prepare static resources and configuration separately, then place them under the directory exposed by Nginx.
The helper scripts in `jbrowse2/scripts/` are intended for this preparation workflow.

## Setup steps
1) Build Hugo static output and sync it to the host path served by Nginx.
2) Prepare JBrowse 2 static resources on the host.
3) Ensure runtime data directories exist for Dash and SequenceServer.
4) Start Dash with `docker compose`.
5) Start SequenceServer with `docker compose`.
6) Include the Nginx snippet and reload Nginx.

## Dash service

The repository already provides:
- `dash/Dockerfile`
- `dash/docker-compose.yml`

Published Docker image:
- `qtulab/miniodp-dash:latest`

Start Dash:

```bash
cd dash
docker compose up -d
```

Or pull the published image directly:

```bash
docker pull qtulab/miniodp-dash:latest
```

Dash expects its converted runtime data under the mounted `dash/data/` directory.
The current compose file mounts:
- source code to `/app`
- `./data` to `/app/data`
- `./logs` to `/app/logs`

## SequenceServer service

The repository uses the official `wurmlab/sequenceserver` image.

Start SequenceServer:

```bash
cd sequenceserver
docker compose up -d
```

The compose file mounts `./data` to `/db` inside the container and starts SequenceServer on port `4567`.

## Nginx hints

Use Nginx on the host to serve static content and reverse proxy the two containerized services.

Repository example:
- `common/nginx/miniodp.conf`

That file routes:
- `/miniodp/` to Hugo static files
- `/miniodp/jbrowse2/` to JBrowse 2 static resources
- `/miniodp/dash/` to `127.0.0.1:3939`
- `/miniodp/sequenceserver/` to `127.0.0.1:4567`

## Environment variables

Dash service variables are defined in `dash/docker-compose.yml`, including:
- `PORT`
- `DASH_HOST`
- `DASH_PREFIX`
- `DASH_DATA_PATH`
- `GUNICORN_*`

Adjust them in the compose file or an external environment file if needed.

## Troubleshooting

- Verify that Hugo output really exists under the directory served by Nginx.
- Verify that JBrowse 2 paths match the Nginx location blocks.
- Verify data mounts and permissions if Dash shows empty studies or failed queries.
- Check Nginx upstream connectivity to `127.0.0.1:3939` and `127.0.0.1:4567`.
- Review container logs when using Docker Compose.

## Demo data

Use separate demo packages for:
- `portal demo`: compact data for Hugo, Dash, JBrowse 2, and SequenceServer
- `pipeline demo`: minimal inputs for `pipeline/`

See these documents for expected structure and usage:

- `docs/PORTAL_GUIDE.md`
- `demo/portal_demo/README.md`
- `demo/pipeline_demo/README.md`
