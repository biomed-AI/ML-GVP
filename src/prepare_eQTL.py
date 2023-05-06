#!/usr/bin/env python3

import argparse, os, sys, time
import numpy as np

from collections import defaultdict
from utils import copen

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-e', '--eqtl-overlap', dest="eqtl", required=True)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    var_eqtls = dict()
    with copen(args.eqtl) as infile:
        for l in infile:
            chrom, eqtl_pos, _, _, start, end, key = l.strip().split('\t')[0:7]
            key = key.split('@@')[0]
            # key = '_'.join([chrom, pos, ref, alt])
            if key not in var_eqtls:
                var_eqtls[key] = defaultdict(int)
            start, end, eqtl_pos = int(start), int(end), int(eqtl_pos) 
            if eqtl_pos < start:
                d = abs(eqtl_pos - start)
            elif eqtl_pos >= end:
                d = abs(eqtl_pos - (end - 1))
            else:
                d = 0

            if d <= 1:
                var_eqtls[key][1] += 1
            elif d <= 10:
                var_eqtls[key][10] += 1
            elif d <= 100:
                var_eqtls[key][100] += 1
            elif d <= 1000:
                var_eqtls[key][1000] += 1
            var_eqtls[key]["weight"] += (1 / np.exp(d/100))

    print('##{}\n##{}'.format(time.asctime(), ' '.join(sys.argv)))
    print('\t'.join(["#key", "eQTL-1", "eQTL-10", "eQTL-100", "eQTL-1000", "eQTL-sum"]))
    with open(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            chrom, pos, key, ref, alt = l.strip().split('\t')[0:5]
            key = key.split('@@')[0]
            # key = '_'.join([chrom, pos, ref, alt])
            if key not in var_eqtls:
                var_eqtls[key] = defaultdict(int)
            print("{}\t{}\t{}\t{}\t{}\t{:.3f}".format(key, var_eqtls[key][1], var_eqtls[key][10], var_eqtls[key][100], var_eqtls[key][1000], var_eqtls[key]["weight"]))


