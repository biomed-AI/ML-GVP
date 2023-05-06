#!/usr/bin/env python3

import argparse, os, sys, time
import numpy as np
import warnings
from utils import copen

# import logging, warnings, json, gzip, pickle
# from collections import defaultdict, OrderedDict
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.stats import pearsonr, spearmanr, ttest_ind, mannwhitneyu
# from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score
def overlap_length(x1, x2, y1, y2):
    """ [x1, x2), [y1, y2) """
    length = 0
    x1, x2, y1, y2 = int(x1), int(x2), int(y1), int(y2)
    if x2 <= y1:
        length = x2 - y1
    elif x1 <= y2:
        length = min(x2, y2) - max(x1, y1)
    else:
        length = y2 - x1
    return length

def distance(x1, x2, y1, y2):
    """ interval distance """
    d = overlap_length(x1, x2, y1, y2)
    return -d

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-isec', required=True)
    return p

CRE_NAME = {
    "dELS": "dELS",
    "dELS,CTCF-bound": "dELS",
    "pELS": "pELS",
    "pELS,CTCF-bound": "pELS",
    "PLS": "PLS",
    "PLS,CTCF-bound": "PLS"
}


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()

    sample_cre = dict()
    with copen(args.isec) as infile:
        for l in infile:
            fields = l.strip().split('\t')
            start1, end1 = int(fields[1]), int(fields[2])
            key = fields[3].split('@@')[0]
            if key not in sample_cre:
                sample_cre[key] = dict()
            start2, end2 = int(fields[5]), int(fields[6])
            cre_type = fields[-1]
            if cre_type not in CRE_NAME:
                warnings.warn(f"unknown cre_type {cre_type}")
                continue
            cre_type = CRE_NAME[cre_type]
            d = max(0, distance(start1, end1, start2, end2))
            assert d >= 0, "{}".format((fields, d, start1, end1, start2, end2))
            d = np.log(1 + d)
            if cre_type not in sample_cre[key]:
                sample_cre[key][cre_type] = d
            else:
                sample_cre[key][cre_type] = min(d, sample_cre[key][cre_type])

    cre_list = sorted(list(set(CRE_NAME.values())))
    print("##{}\n##{}\n#key\t{}".format(time.asctime(), ' '.join(sys.argv), '\t'.join(cre_list)))
    with open(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('@@')[0]
            if key in sample_cre:
                print("{}\t{}".format(key, '\t'.join([
                    '{:.3g}'.format(sample_cre[key][r]) if r in sample_cre[key] else '10' for r in cre_list
                    ])))
            else:
                warnings.warn("missing for {}".format(key))
                print("{}\t{}".format(key, '\t'.join(['10' for r in cre_list])))

            # print("{}\t{}\t{}".format(key, len(var_rbp_count["K562"][key]), len(var_rbp_count["HepG2"][key])))


