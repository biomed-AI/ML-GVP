#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
import pandas as pd
from utils import copen, run_bash

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-r', '--ref', required=True, help="reference genome fasta")
    p.add_argument("--overwrite", action="store_true")
    return p


if __name__ == "__main__":
    args = get_args().parse_args()

    fai = args.ref + ".fai"

    if not os.path.exists(fai):
        run_bash(f"samtools faidx {args.ref}")
    
    header = list()
    header.append('##fileformat=VCFv4.0')
    header.append('##FILTER=<ID=PASS,Description="All filters passed">')
    fasta = os.path.basename(args.ref)
    with open(fai) as infile:
        for l in infile:
            chrom, size = l.strip().split('\t')[:2]
            header.append(f'##contig=<ID={chrom},length={size},assembly={fasta}>')
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    cache = list()
    if args.overwrite:
        cache.append('\n'.join(header))
    else:
        print('\n'.join(header))

    with copen(args.vcf) as infile:
        for l in infile:
            if l.startswith("#"):
                continue
            else:
                if args.overwrite:
                    cache.append(l.strip())
                else:
                    print(l, end='')
    if args.overwrite:
        with copen(args.vcf, mode='wt') as out:
            out.write("\n".join(cache))
