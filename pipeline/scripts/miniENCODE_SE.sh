#!/bin/bash
set -euo pipefail

# miniENCODE preprocessing pipeline
# Script Purpose: Identify Super Enhancers using ROSE.

# Set the path to ROSE
dir=$PWD
PATHTO="/usr/local/share/ROSE-1.3.2"
PYTHONPATH="$PATHTO/lib"
export PYTHONPATH
export PATH="$PATH:$PATHTO/bin"
ROSE="$PATHTO/bin/ROSE_main.py"

# Set the number of threads to use
threads="${THREADS:-2}"

# Define the working directories
inputDir="${dir}/data/alignment/"
outDir="${dir}/data/SE/"

# Change to the output directory
cd "$outDir" || exit
chmod +x "$ROSE"
mkdir -p "${outDir}/tmp"

# Define the paths to BAM files and narrowPeak file for the current sample
sp="${SE_SAMPLE:-Heart}"
bam="${SE_SIGNAL_BAM:-${inputDir}/ch_H3K27ac_heart_2020nat_chr1.final.bam}"
inputbam="${SE_INPUT_BAM:-${inputDir}/ch_input_heart_2020nat_chr1.final.bam}"
narrowPeakFile="${SE_PEAK_FILE:-${dir}/data/peak/ch_H3K27ac_heart_2020nat_chr1_peaks.narrowPeak}"
refseq="${SE_REFSEQ:-${dir}/data/reference/GRCz11_ucsc.refseq}"

[ -f "${bam}" ] || { echo "Missing signal BAM: ${bam}" >&2; exit 1; }
[ -f "${inputbam}" ] || { echo "Missing input BAM: ${inputbam}" >&2; exit 1; }
[ -f "${narrowPeakFile}" ] || { echo "Missing narrowPeak file: ${narrowPeakFile}" >&2; exit 1; }
[ -f "${refseq}" ] || { echo "Missing refseq file: ${refseq}" >&2; exit 1; }

# Reheader BAM files to replace chromosome names
samtools view -@ "$threads" -H "$bam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
    | samtools reheader - "$bam" > "${outDir}/tmp/ch_H3K27ac_${sp}_chr.bam"
samtools index -@ "$threads" "${outDir}/tmp/ch_H3K27ac_${sp}_chr.bam"

samtools view -@ "$threads" -H "$inputbam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
    | samtools reheader - "$inputbam" > "${outDir}/tmp/ch_H3K27ac_${sp}_input_chr.bam"
samtools index -@ "$threads" "${outDir}/tmp/ch_H3K27ac_${sp}_input_chr.bam"

# Modify the narrowPeak file to have "chr" prefixes
sed 's/^/chr/g' "$narrowPeakFile" > "${outDir}/tmp/ch_H3K27ac_${sp}_narrowPeak.bed"

# Run ROSE with the specified parameters
mkdir -p "${outDir}/${sp}"
python "$ROSE" --custom="$refseq" -i "${outDir}/tmp/ch_H3K27ac_${sp}_narrowPeak.bed" \
    -r "${outDir}/tmp/ch_H3K27ac_${sp}_chr.bam" \
    -c "${outDir}/tmp/ch_H3K27ac_${sp}_input_chr.bam" \
    -o "${outDir}/${sp}/" \
    -s 12500 -t 2500 2>"${outDir}/${sp}/${sp}.log"

# Remove temporary directory
rm -rf "${outDir}/tmp"
