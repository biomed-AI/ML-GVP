#!/usr/bin/env python3

import argparse
import os
from collections import OrderedDict
import sys
import numpy as np
import pandas as pd
from utils import copen


def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("vcf")
    # p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    args = get_args().parse_args()

    std_vcf = list()

    with copen(args.vcf) as infile:
        idx = 0
        for l in infile:
            if l.startswith("#"):
                continue
            chrom, position, name, ref, alt = l.strip().split('\t')[0:5]
            key = f"{chrom}_{position}_{ref}_{alt}"
            std_vcf.append((
                chrom, 
                int(position), 
                f"{idx}@{key}@@{name}",
                ref, alt
            ))
            idx += 1

    std_vcf = sorted(std_vcf, key=lambda x:(x[0], x[1]))

    for chrom, pos, name, ref, alt in std_vcf:
        print(f"{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t.\t.\t.")
