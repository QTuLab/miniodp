# Pipeline Demo

This directory documents the separate `pipeline_demo` package used to validate the analysis workflow in `pipeline/`.

## Scope

The demo package is organized around zebrafish chromosome 1 and covers:

- download and QC
- trimming
- reference preparation
- RNA-seq
- ATAC-seq
- ChIP-seq
- BS-seq
- ELS
- super-enhancer prediction
- GRN inference

## Download

- file: `pipeline_demo.tar.gz`
- url: `https://tulab.genetics.ac.cn/~qtu/miniodp/pipeline_demo.tar.gz`
- size: `8.4G`
- SHA256: `de8db26fe01f6df91c27d14665438397a48f4ff4b2872b69977649e8c66488e6`

The checksum list is published at:

- `https://tulab.genetics.ac.cn/~qtu/miniodp/SHA256SUMS.txt`

## Layout

```text
pipeline_demo/
└── data/
    ├── info/
    ├── reference/
    ├── srafile/
    ├── fastq/
    ├── qc/
    ├── alignment/
    ├── quant/
    ├── peak/
    ├── ELS/
    ├── SE/
    └── GRN/
```

This package is released as a validation bundle rather than a minimal input bundle.
It intentionally keeps representative intermediate outputs so that each workflow stage can be checked without rerunning the full pipeline.

## Core inputs

The package should contain:

- `data/info/`
  - `id_download.txt`
  - `id_rna.txt`
  - `id_atac.txt`
  - `id_chip.txt`
  - `id_bs.txt`
  - `namelist_rna.csv`
  - `namelist_atac.csv`
  - `namelist_chip.csv`
  - `namelist_chip_callpeak.csv`
  - `namelist_bs.csv`
- `data/reference/`
  - `GRCz11.fa`
  - `GRCz11_chr1_dna.fa`
  - `GRCz11_chr1_100.gtf`
  - `GRCz11_chr1_100.gff3`
  - `GRCz11_chr1_100.bed`
  - `GRCz11_chr1.chrom.sizes`
  - `GRCz11_ucsc.refseq`
  - `GRCz11.annotation.gene.ensembl.bed`
  - `GRCz11_TSS_Basic.bed`
  - `GRCz11_TSS_Basic4K.bed`
  - `GRCz11_TFDB/`
- analysis directories:
  - `data/srafile/`
  - `data/fastq/`
  - `data/qc/`
  - `data/alignment/`
  - `data/quant/`
  - `data/peak/`
  - `data/ELS/`
  - `data/SE/`
  - `data/GRN/`

## Reference sample groups

Core sample groups used in the compact demo:

- RNA:
  - `tissue_heart_2020nat_chr1`
  - `st_24hpf_2011cell_chr1`
- ATAC:
  - `at_brain_2020nat_chr1`
- ChIP:
  - `ch_H3K27ac_heart_2020nat_chr1`
  - `ch_input_heart_2020nat_chr1`
- BS:
  - `bs_pgc_2019nc_chr1`
- ELS:
  - `Blastula`
  - `Blood`
- GRN:
  - `Blood`

## Validation path

To validate the current repository with this package:

1. Unpack the package so that `pipeline_demo/data/` becomes `pipeline/data/`.
2. Build the runtime image from `pipeline/Dockerfile`.
3. Start the runtime shell with `cd pipeline && docker compose run --rm pipeline bash`.
4. Inside the container, run `cd /work/pipeline`.
5. Run the scripts in the order documented in [pipeline/README.md](../../pipeline/README.md).

## Compatibility notes

- The demo package is derived from legacy miniENCODE demonstration materials and was reorganized for the current repository layout.
- The current `miniENCODE_GRN.sh` script matches the legacy demo filenames under `data/GRN/input/`.
- GRN inference uses the whole-genome reference file `data/reference/GRCz11.fa`.
- Full end-to-end validation depends on the container environment defined by `pipeline/Dockerfile`.
