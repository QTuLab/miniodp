# Pipeline Guide

This document focuses on the analysis workflow in `pipeline/`.
For the full project overview, see the repository root [README.md](../README.md).

## Scope

The current pipeline scripts cover:

- download and QC
- trimming
- reference preparation
- RNA-seq
- ATAC-seq
- ChIP-seq
- BS-seq
- enhancer-like signature analysis
- super-enhancer prediction
- gene regulatory network inference

The repository provides a containerized runtime for these steps through:

- [Dockerfile](Dockerfile)
- [docker-compose.yml](docker-compose.yml)

Published Docker image:

- `qtulab/miniodp-pipeline:latest`

## Runtime model

Build the image:

```bash
cd pipeline
docker build -t qtulab/miniodp-pipeline:latest .
```

Or pull the published image directly:

```bash
docker pull qtulab/miniodp-pipeline:latest
```

Open an interactive shell:

```bash
cd pipeline
docker compose run --rm pipeline bash
```

The compose file mounts the parent directory of `pipeline/` into `/work` inside the container. If needed, override the mount root:

```bash
cd pipeline
PIPELINE_WORKDIR=/absolute/path/to/project docker compose run --rm pipeline bash
```

## Expected directory layout

The workflow expects its working data under `pipeline/data/`.

```text
pipeline/
├── data/
│   ├── info/
│   ├── reference/
│   ├── srafile/
│   ├── fastq/
│   ├── qc/
│   ├── alignment/
│   ├── quant/
│   ├── peak/
│   ├── ELS/
│   ├── SE/
│   └── GRN/
└── scripts/
```

See the separate demo package documentation for the expected contents of `pipeline/data/`:

- [../demo/pipeline_demo/README.md](../demo/pipeline_demo/README.md)

Current demo download:

- `https://tulab.genetics.ac.cn/~qtu/miniodp/pipeline_demo.tar.gz`
- SHA256: `de8db26fe01f6df91c27d14665438397a48f4ff4b2872b69977649e8c66488e6`

## Workflow order

Run all wrapper scripts from inside `pipeline/`.
If you start an interactive container shell from the repository root, change into `/work/pipeline` before running the commands below.

### 1. Download and QC

Purpose:
- fetch public sequence data
- run initial QC

Entry point:

```bash
bash scripts/miniENCODE_pre_download_qc.sh
```

Main inputs:
- `pipeline/data/info/id_download.txt`

### 2. Trimming

Purpose:
- trim reads when required by QC

Entry point:

```bash
bash scripts/miniENCODE_pre_trim_qc.sh
```

### 3. Reference preparation

Purpose:
- prepare alignment indexes and reference files

Entry point:

```bash
bash scripts/miniENCODE_pre_makeRef.sh
```

Typical reference inputs:
- `pipeline/data/reference/GRCz11_chr1_dna.fa`
- `pipeline/data/reference/GRCz11_chr1_100.gtf`

### 4. RNA-seq

Purpose:
- alignment
- merge
- quantification

Entry point:

```bash
bash scripts/miniENCODE_pre_RNA.sh
```

Main control files:
- `pipeline/data/info/id_rna.txt`
- `pipeline/data/info/namelist_rna.csv`

### 5. ATAC-seq

Purpose:
- alignment
- peak calling

Entry point:

```bash
bash scripts/miniENCODE_pre_ATAC.sh
```

Main control files:
- `pipeline/data/info/id_atac.txt`
- `pipeline/data/info/namelist_atac.csv`

Peak-calling note:
- The current script calls MACS2 with `-f BAMPE` for paired-end ATAC-seq and keeps the demo genome size `-g 1.1e9`.
- If your data are single-end or use a different genome, update the script or move these settings into a config file before running the workflow.

### 6. ChIP-seq

Purpose:
- alignment
- peak calling

Entry point:

```bash
bash scripts/miniENCODE_pre_ChIP.sh
```

Main control files:
- `pipeline/data/info/id_chip.txt`
- `pipeline/data/info/namelist_chip.csv`
- `pipeline/data/info/namelist_chip_callpeak.csv`

### 7. BS-seq

Purpose:
- alignment and methylation-oriented processing

Entry point:

```bash
bash scripts/miniENCODE_pre_BS.sh
```

Main control files:
- `pipeline/data/info/id_bs.txt`
- `pipeline/data/info/namelist_bs.csv`

### 8. Enhancer-like signatures

Purpose:
- derive rDAR regions
- calculate z-scores
- classify ELS groups

Entry points:

```bash
bash scripts/miniENCODE_ELS_rDAR.sh
bash scripts/miniENCODE_ELS_zscore.sh
bash scripts/miniENCODE_ELS_classifyELS.sh
```

Typical inputs:
- `pipeline/data/reference/GRCz11_TSS_Basic.bed`
- `pipeline/data/reference/GRCz11_TSS_Basic4K.bed`
- `pipeline/data/ELS/`

### 9. Super-enhancer prediction

Purpose:
- run ROSE-based super-enhancer calling from H3K27ac and control alignments

Entry point:

```bash
bash scripts/miniENCODE_SE.sh
```

Typical inputs:
- `pipeline/data/reference/GRCz11_ucsc.refseq`
- `pipeline/data/alignment/ch_H3K27ac_heart_2020nat_chr1.final.bam`
- `pipeline/data/alignment/ch_input_heart_2020nat_chr1.final.bam`

### 10. GRN inference

Purpose:
- infer enhancer-gene regulatory networks from ATAC, H3K27ac, ELS, expression, and TF motif resources

Entry point:

```bash
bash scripts/miniENCODE_GRN.sh
```

Typical inputs:
- `pipeline/data/GRN/input/at_Blood.bam`
- `pipeline/data/GRN/input/ch_H3K27ac_Blood.bam`
- `pipeline/data/GRN/input/els_Blood.bed`
- `pipeline/data/GRN/input/tpm_Blood.txt`
- `pipeline/data/reference/GRCz11_TFDB/`

## Demo validation path

To validate the current repository with the separate pipeline demo package:

1. Unpack the demo so that its `data/` directory becomes `pipeline/data/`.
2. Build the image from `pipeline/Dockerfile`, or pull `qtulab/miniodp-pipeline:latest` if you only need the published runtime.
3. Start the container with `docker compose run --rm pipeline bash`.
4. Inside the container, run `cd /work/pipeline`.
5. Run the scripts in workflow order from `/work/pipeline`.

## BulkMulti conversion notes

If you convert pipeline outputs into BulkMulti format with `dash/scripts/bulkmulti_convert_from_pipeline.py`, the default path inference expects:

- `data/GRN/<sample_id>/<sample_id>.network.txt`
- `data/ELS/<sample_id>/ELS_<sample_id>.bed`
- `data/GRN/input/tpm_<sample_id>.txt`
- `data/alignment/<sample_id>.bam.bw`

You can override any of these paths in the metadata CSV. The converter log will show whether each path came from the default rule or a metadata override.

The demo package is intended to validate path resolution, software availability, and script compatibility with a compact zebrafish chromosome 1 dataset.

## Notes on compatibility

- The public demo package is derived from legacy miniENCODE materials and was reorganized for the current repository layout.
- Full end-to-end verification depends on the container environment provided by `pipeline/Dockerfile`.
- The current image includes the validated toolchain needed by the scripts, including HISAT2, Bowtie 2, StringTie, SAMtools, BEDTools, FastQC, Trim Galore, Picard, SRA Toolkit, Bismark, deepTools, MACS2, ROSE, Graphviz, and ANANSE.
