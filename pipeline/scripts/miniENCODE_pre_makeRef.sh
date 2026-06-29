#!/bin/bash
set -euo pipefail

# miniENCODE preprocessing pipeline
# Script Purpose: Build HISAT2, Bowtie2, and Bismark indexes for genomic reference data.

# Define common variables
dir=$PWD
scriptDir="${dir}/scripts"
refDir="${dir}/data/reference/"
fa="${REFERENCE_FASTA:-${refDir}/GRCz11_chr1_dna.fa}"
gtf="${REFERENCE_GTF:-${refDir}/GRCz11_chr1_100.gtf}"
thread="${THREADS:-4}"

[ -f "${fa}" ] || { echo "Missing reference FASTA: ${fa}" >&2; exit 1; }
[ -f "${gtf}" ] || { echo "Missing reference GTF: ${gtf}" >&2; exit 1; }

# Build HISAT2 index
indexDir="${refDir}/index_hisat2/"
rm -rf "${indexDir}"
mkdir -p "${indexDir}"
cd "${indexDir}" || exit

python "${scriptDir}/hisat2_extract_splice_sites.py" "${gtf}" > "${indexDir}/GRCz11.ss"
python "${scriptDir}/hisat2_extract_exons.py" "${gtf}" > "${indexDir}/GRCz11.exon"
hisat2-build -p "${thread}" --ss "${indexDir}/GRCz11.ss" --exon "${indexDir}/GRCz11.exon" "${fa}" GRCz11

# Build Bowtie2 index
indexDir="${refDir}/index_bowtie2/"
rm -rf "${indexDir}"
mkdir -p "${indexDir}"
cd "${indexDir}" || exit

bowtie2-build "${fa}" GRCz11

# Build Bismark index
indexDir="${refDir}/index_bismark/"
rm -rf "${indexDir}"
mkdir -p "${indexDir}"
cd "${indexDir}" || exit
bisDir=$(dirname "$(which bowtie2)")
ln -s "${fa}" "${indexDir}"

bismark_genome_preparation --path_to_aligner "${bisDir}" --verbose "${indexDir}"
