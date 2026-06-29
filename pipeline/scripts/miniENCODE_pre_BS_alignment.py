#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script aligns BS-seq data using Bismark.

import sys
import getopt
import logging
import miniENCODE_function as mf

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def usage():
    log.error("Error!")
    log.error("Usage: " + sys.argv[0] + " -i|--input <inputSRAListfile> -o|--output <outpathway> -t|--threads <numberOfThreads>")
    sys.exit(2)

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:t:', ['input=', 'output=', 'threads='])
except getopt.GetoptError:
    usage()

inputfile = ''
outpathway = ''
threads = 1

for opt, arg in opts:
    if opt in ('-i', '--input'):
        inputfile = arg
    elif opt in ('-o', '--output'):
        outpathway = arg
    elif opt in ('-t', '--threads'):
        threads = arg

if not inputfile or not outpathway:
    usage()

mf.require_file(inputfile, 'input SRA list file')

def bs_seq_gz(dataDir, bisindex, n, sra):
    fq_dir = f'{dataDir}/fastq/'
    align_dir = f'{dataDir}/alignment/'
    layout = sra[1]
    sra = sra[0]
    n = str(n)
    if layout == 'paired':
        fq1 = f'{fq_dir}{sra}_1_val_1.fq.gz'
        fq2 = f'{fq_dir}{sra}_2_val_2.fq.gz'
        mf.run_shell(f'bismark --parallel {n} -o {align_dir} --genome {bisindex} -1 {fq1} -2 {fq2}')
        mf.run_shell(f'deduplicate_bismark -o {sra} --output_dir {align_dir} --bam {align_dir}{sra}_1_val_1_bismark_bt2_pe.bam')
    elif layout == 'single':
        mf.run_shell(f'bismark --parallel {n} -o {align_dir} --genome {bisindex} {fq_dir}{sra}_trimmed.fq.gz')
        mf.run_shell(f'deduplicate_bismark -o {sra} --output_dir {align_dir} --bam {align_dir}{sra}_trimmed_bismark_bt2.bam')

mf.run_shell(f'mkdir -p {outpathway}/alignment')

bisindex = f'{outpathway}/reference/index_bismark/'
mf.require_dir(bisindex, 'Bismark genome index directory')

with open(inputfile, 'r') as SraListFile:
    for line in SraListFile:
        line = line.strip()
        log.info(line)
        line = line.split(',')
        bs_seq_gz(outpathway, bisindex, threads, line)
