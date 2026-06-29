# Dash Application

## Purpose
This Dash app delivers interactive multi-omics views (gene info, bulk RNA, scRNA, scATAC, BulkMulti) with preprocessed data for fast queries.

## Structure
- `app.py`: entry point
- `assets/`: static assets and shared styles
- `backend/`: API routes and data accessors
- `internal/`: internal utilities shared across routes
- `scripts/`: data preparation helpers
- `data/`: converted datasets (ignored in git; populated externally)

## Local development
1) Environment: `mamba env create -f environment.yml && mamba activate miniodp-dash`
2) Install runtime deps only: `pip install -r requirements.txt`
3) Run: `python app.py`

## Container image

Published Docker image:

- `qtulab/miniodp-dash:latest`

Pull it with:

```bash
docker pull qtulab/miniodp-dash:latest
```

## Notes
- Portal integration and deployment notes are documented in `../docs/PORTAL_GUIDE.md` and `../docs/DEPLOYMENT_GUIDE.md`.
- Data conversion and new-species setup are documented in `../docs/DATA_PREPARATION_GUIDE.md` and `../docs/ONBOARDING_GUIDE.md`.
