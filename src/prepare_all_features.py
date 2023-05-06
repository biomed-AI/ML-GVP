#!/usr/bin/env python3

import argparse, os, sys, time
from tqdm import tqdm
from functools import partial
print = partial(print, flush=True)

# import logging, warnings, json, gzip, pickle
# from collections import defaultdict, OrderedDict
import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.stats import pearsonr, spearmanr, ttest_ind, mannwhitneyu
# from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('inputs', nargs='+')
    p.add_argument('-o', required=True)

    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    df = None
    pbar = tqdm(args.inputs)
    for fn in pbar:
        pbar.set_postfix_str(f"reading {fn}", refresh=False)
        print("{}".format(fn))
        if df is None:
            df = pd.read_csv(fn, delimiter='\t', index_col=0, skiprows=2, dtype=str)
        else:
            df = pd.concat((df, pd.read_csv(fn, delimiter='\t', index_col=0, skiprows=2, dtype=str)), axis=1)


    df.to_csv(args.o, sep='\t', na_rep="NaN")
