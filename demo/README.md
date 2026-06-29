# miniodp Demonstration Data Package

This dataset contains two separate demonstration archives for the `miniodp`
framework. They serve different purposes and are intended to be downloaded
independently.

## Files

### 1. `portal_demo.tar.gz`

A compact portal-side example derived from a zebrafish deployment.

It contains a small but working subset for:

- Hugo species display configuration
- Dash runtime data
- JBrowse 2 species bundle
- SequenceServer example database

Use this archive to validate portal deployment and cross-component wiring.

SHA-256:

`85f5d2502744bc740602aeaecc3e0d9b5665c080fa5bab923d549107377d223c`

### 2. `pipeline_demo.tar.gz`

A compact pipeline-side example centered on zebrafish chromosome 1.

It contains:

- reference files
- control tables
- FASTQ and intermediate outputs
- representative inputs for downstream integrative analyses

Use this archive to validate the analysis workflow and runtime environment.

SHA-256:

`de8db26fe01f6df91c27d14665438397a48f4ff4b2872b69977649e8c66488e6`

## Notes

- These files are demonstration packages, not full production datasets.
- The two archives are intentionally distributed separately so users can
  download only the part they need.
- The corresponding source code and detailed documentation are available in the
  public repository: `https://github.com/qtulab/miniodp`

