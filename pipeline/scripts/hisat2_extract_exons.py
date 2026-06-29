#!/usr/bin/env python3

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#

from argparse import ArgumentParser, FileType
from collections import defaultdict as dd
from sys import exit


def extract_exons(gtf_file, verbose=False):
    genes = dd(list)
    trans = {}

    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "#" in line:
            line = line.split("#")[0].strip()

        try:
            chrom, source, feature, left, right, score, strand, frame, values = (
                line.split("\t")
            )
        except ValueError:
            continue
        left, right = int(left), int(right)

        if feature != "exon" or left >= right:
            continue

        values_dict = {}
        for attr in values.split(";")[:-1]:
            attr, _, val = attr.strip().partition(" ")
            values_dict[attr] = val.strip('"')

        if "gene_id" not in values_dict or "transcript_id" not in values_dict:
            continue

        transcript_id = values_dict["transcript_id"]
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
            genes[values_dict["gene_id"]].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    for tran, (chrom, strand, exons) in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]

    tmp_exons = set()
    for chrom, strand, texons in trans.values():
        for exon in texons:
            tmp_exons.add((chrom, exon[0], exon[1], strand))
    tmp_exons = sorted(tmp_exons)
    if not tmp_exons:
        return

    exons = [tmp_exons[0]]
    for exon in tmp_exons[1:]:
        prev_exon = exons[-1]
        if exon[0] != prev_exon[0]:
            exons.append(exon)
            continue
        assert prev_exon[1] <= exon[1]
        if prev_exon[2] < exon[1]:
            exons.append(exon)
            continue

        if prev_exon[2] < exon[2]:
            strand = prev_exon[3]
            if strand not in "+-":
                strand = exon[3]
            exons[-1] = (prev_exon[0], prev_exon[1], exon[2], strand)

    for chrom, left, right, strand in exons:
        print(f"{chrom}\t{left - 1}\t{right - 1}\t{strand}")

    if verbose:
        return


if __name__ == "__main__":
    parser = ArgumentParser(description="Extract exons from a GTF file")
    parser.add_argument(
        "gtf_file",
        nargs="?",
        type=FileType("r"),
        help='input GTF file (use "-" for stdin)',
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="also print some statistics to stderr",
    )

    args = parser.parse_args()
    if not args.gtf_file:
        parser.print_help()
        exit(1)
    extract_exons(args.gtf_file, args.verbose)
