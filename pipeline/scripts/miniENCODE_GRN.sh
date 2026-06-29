#!/bin/bash
set -euo pipefail

# miniENCODE preprocessing pipeline
# Script Purpose: Perform gene regulatory network (GRN) analysis using ANANSE.

# Define reference directory and input/output directories
dir=$PWD
refDir="${dir}/data/reference/"
inputDir="${dir}/data/GRN/input/"
outDir="${dir}/data/GRN/"
reftfbs="${GRN_TF_MOTIF:-${refDir}/GRCz11_TFDB/danRer11.gimme.vertebrate.v5.0.pfm}"
reference_fasta="${GRN_REFERENCE_FASTA:-${refDir}/GRCz11.fa}"
reference_annotation="${GRN_REFERENCE_ANNOTATION:-${refDir}/GRCz11.annotation.gene.ensembl.bed}"

# Set the number of threads
thread="${THREADS:-1}"
samples="${GRN_SAMPLES:-Blood}"

for sp in ${samples}
do
    mkdir -p "${outDir}/${sp}"

    # Define input files for ANANSE
    atacbam="${inputDir}/at_${sp}.bam"
    chipbam="${inputDir}/ch_H3K27ac_${sp}.bam"
    elsfile="${inputDir}/els_${sp}.bed"
    tpmfile="${inputDir}/tpm_${sp}.txt"

    [ -f "${atacbam}" ] || { echo "Missing ATAC BAM: ${atacbam}" >&2; exit 1; }
    [ -f "${chipbam}" ] || { echo "Missing H3K27ac BAM: ${chipbam}" >&2; exit 1; }
    [ -f "${elsfile}" ] || { echo "Missing ELS BED: ${elsfile}" >&2; exit 1; }
    [ -f "${tpmfile}" ] || { echo "Missing TPM file: ${tpmfile}" >&2; exit 1; }
    [ -f "${reference_fasta}" ] || { echo "Missing reference FASTA: ${reference_fasta}" >&2; exit 1; }
    [ -f "${reference_annotation}" ] || { echo "Missing reference annotation: ${reference_annotation}" >&2; exit 1; }
    [ -f "${reftfbs}" ] || { echo "Missing TF motif file: ${reftfbs}" >&2; exit 1; }

    # Run ANANSE binding analysis
    ananse binding -A "$atacbam" \
                   -H "$chipbam" \
                   -r "$elsfile" \
                   -g "${reference_fasta}" \
                   -p "$reftfbs" \
                   -o "${outDir}/${sp}" \
                   -n "$thread"

    # Convert ANANSE binding output to TSV format
    ananse view "${outDir}/${sp}/binding.h5" -o "${outDir}/${sp}/binding.tsv"

    # Run ANANSE network analysis
    ananse network -g "${reference_fasta}" \
                   -a "${reference_annotation}" \
                   -e "${tpmfile}" \
                   -o "${outDir}/${sp}/${sp}.network.txt" \
                   -n "$thread" \
                   "${outDir}/${sp}/binding.h5" > "${outDir}/${sp}/log" 2>&1 &
done

wait
