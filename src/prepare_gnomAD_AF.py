#!/usr/bin/env python3

import argparse, os, sys, time
import numpy as np
from typing import Any, Dict, List, Union
from functools import partial
print = partial(print, flush=True)

# import logging, warnings, json, gzip, pickle
# from collections import defaultdict, OrderedDict
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.stats import pearsonr, spearmanr, ttest_ind, mannwhitneyu
# from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('--af', required=True, nargs='+')

    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)

    var_af = dict()
    for fn in args.af:
        with open(fn) as infile:
            for l in infile:
                if not l.startswith("chr"):
                    continue
                chrom, hg38_pos1, _, ref1, alt1, _, _, af, _, hg38_pos2, name, ref2, alt2, _, _, _, _ = l.strip().split('\t')
                if hg38_pos1 != hg38_pos2:
                    continue
                if ref1 != ref2 or alt1 != alt2:
                    continue
                key = name.split('@@')[0]
                var_af[key] = af.split(',')[0]
                if var_af[key] == '.':
                    var_af[key] = "NaN"

    print("##{}\n##{}\n#key\tgnomAD.v3.AF".format(time.asctime(), ' '.join(sys.argv)))
    with open(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            name = l.split('\t')[2].split('@@')[0]
            if name in var_af:
                print("{}\t{}".format(name, var_af[name]))
            else:
                print("{}\t0.0000".format(name))


