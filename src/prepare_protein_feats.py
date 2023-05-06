#!/usr/bin/env python3

import argparse, os, sys, time
import numpy as np
from typing import Any, Dict, List, Union
from functools import partial
print = partial(print, flush=True)

import logging, warnings, json, gzip, pickle
# from collections import defaultdict, OrderedDict
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.stats import pearsonr, spearmanr, ttest_ind, mannwhitneyu
# from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-p', '--pssm-disorder', required=True, help="preprocessed dict")
    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)

    d = pickle.load(gzip.open(args.pssm_disorder, 'rb'))
    

    print("##{}\n##{}\n#key\t{}\t{}\t{}".format(
        time.asctime(), 
        ' '.join(sys.argv), 
        '\t'.join([
            "disorder-max", "disorder-mean", 
            "disorder.3aa-max", "disorder.3aa-mean", 
            "disorder.10aa-max", "disorder.10aa-mean"
        ]),
        '\t'.join(["PSSM-{}".format(aa) for aa in "ARNDCQEGHILKMFPSTWYV"]),
        "PSSM-change.sumlog"
    ))
    with open(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('|')[0]
            name = key.split("@@")[0]
            if key not in d:
                print("{}\t{}".format(name, '\t'.join(["NA" for _ in range(27)])))
            else:
                try:
                    print("{}\t{}\t{}\t{:.3g}".format(
                        name, 
                        '\t'.join([
                            "{:.3g}".format(np.max(d[key]["disorder"][0])),
                            "{:.3g}".format(np.mean(d[key]["disorder"][0])),
                            "{:.3g}".format(np.max(d[key]["disorder"][1])),
                            "{:.3g}".format(np.mean(d[key]["disorder"][1])),
                            "{:.3g}".format(np.max(d[key]["disorder"][2])),
                            "{:.3g}".format(np.mean(d[key]["disorder"][2])),
                        ]),
                        '\t'.join(["{:.1f}".format(x) for x in d[key]["pssm"].mean(axis=0)]),
                        np.log(1 + np.sum(np.abs(d[key]["pssm_change"])))
                    ))
                except:
                    raise RuntimeError("{}".format((key, d[key])))

