#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script merges ChIP-seq data and calls peaks using MACS2.

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
    log.error("Usage: " + sys.argv[0] + " -s|--sample <sampleListFile> -c|--callpeak <callpeakFile> -o|--output <outpathway> -t|--threads <numberOfThreads>")
    sys.exit(2)

try:
    opts, args = getopt.getopt(sys.argv[1:], 's:c:o:t:', ['sample=', 'callpeak=', 'output=', 'threads='])
except getopt.GetoptError:
    usage()

samplefile = ''
callpeakfile = ''
outpathway = ''
threads = 1

for opt, arg in opts:
    if opt in ('-s', '--sample'):
        samplefile = arg
    elif opt in ('-c', '--callpeak'):
        callpeakfile = arg
    elif opt in ('-o', '--output'):
        outpathway = arg
    elif opt in ('-t', '--threads'):
        threads = arg

if not samplefile or not callpeakfile or not outpathway:
    usage()

mf.require_file(samplefile, 'sample list file')
mf.require_file(callpeakfile, 'callpeak list file')

def chip_seq_merge(dataDir, relist, n):
    align_dir = f'{dataDir}/alignment/'
    bam_file = f'{align_dir}{relist[0]}.final.bam'
    n = str(n)
    if len(relist) == 2:
        log.info('Rename data')
        mf.run_shell(f'mv {align_dir}{relist[1]}.final.bam {bam_file}')
    elif len(relist) > 2:
        nfile = len(relist)
        log.info(f'Merge {nfile} files')
        merge_cmd = ' '.join([f'{align_dir}{file}.final.bam' for file in relist])
        mf.run_shell(f'samtools merge -f -@ {n} {merge_cmd}')
    mf.run_shell(f'samtools sort -@ {n} -o {bam_file} {bam_file}')
    mf.run_shell(f'samtools index -@ {n} {bam_file}')
    mf.run_shell(f'bamCoverage -bs 1 --normalizeUsing RPKM --numberOfProcessors {n} -b {bam_file} -o {bam_file}.bw')

def chip_callpeak(dataDir, relist, n):
    align_dir = f'{dataDir}/alignment/'
    peak_dir = f'{dataDir}/peak/'
    log.info('Peak calling')
    if len(relist) == 1:
        mf.run_shell(f'macs2 callpeak --outdir {peak_dir} --keep-dup all -g 1.1e9 -t {align_dir}{relist[0]}.final.bam -n {relist[0]}')
    elif len(relist) == 2:
        mf.run_shell(f'macs2 callpeak --outdir {peak_dir} --keep-dup all -g 1.1e9 -t {align_dir}{relist[0]}.final.bam -c {align_dir}{relist[1]}.final.bam -n {relist[0]}')

with open(samplefile, 'r') as SampleListFile:
    for line in SampleListFile:
        line = line.strip()
        log.info(line)
        relist = line.split(',')
        chip_seq_merge(outpathway, relist, threads)

with open(callpeakfile, 'r') as CallPeakFile:
    for line in CallPeakFile:
        line = line.strip()
        log.info(line)
        relist = line.split(',')
        chip_callpeak(outpathway, relist, threads)
