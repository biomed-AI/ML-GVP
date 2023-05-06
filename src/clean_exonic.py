#!/usr/bin/env python3

import argparse, os, sys, time
import logging, warnings, json, gzip, pickle


#from collections import defaultdict, OrderedDict

import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#from scipy.stats import pearsonr, spearmanr
from typing import Any, Dict, List, Union
from functools import partial
print = partial(print, flush=True)

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # p.add_argument('-vcf', required=True)
    p.add_argument('-exonic', required=True)
    p.add_argument('-m', '--tx2uniprot', required=True, help="transcript id to uniprot id mapping, preferring txs with uniprot sequencs")

    #p.add_argument('--seed', type=int, default=2020)
    return p

def remove_version(name):
    if name.endswith("PAR_Y"):
        suffix = "_PAR_Y"
    else:
        suffix = ""
    return name.split('.')[0] + suffix


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)

    # line2tx = list()
    # with open(args.vcf) as infile:
    #     for l in infile:
    #         try:
    #             tx_id = l.split('\t')[2].split('|')[-7]
    #         except:
    #             raise RuntimeError("{}".format(l))
    #         assert tx_id.startswith('NM_') or tx_id.startswith("NR_"), tx_id
    #         line2tx.append(tx_id)
    #         del tx_id
    tx2uniprot_raw = json.load(open(args.tx2uniprot))
    tx2uniprot = dict()
    for tx_id in tx2uniprot_raw:
        uniprot_id, score, prot_id, gene_id, gene_name = tx2uniprot_raw[tx_id]
        tx2uniprot[remove_version(tx_id)] = score

    
    with open(args.exonic) as infile:
        for l in infile:
            fields = l.strip().split('\t')
            lid = int(fields[0].replace("line", '')) - 1
            exonic_info = fields[2].rstrip(',').split(',')
            einfo = None
            success = False
            # if len(exonic_info) > 0:
            changes = list()
            for einfo in exonic_info:
                try:
                    changes.append((einfo, tx2uniprot[einfo.split(':')[1]]))
                    success = True
                except:
                    # warnings.warn("Missing value for {} in {}".format(einfo, exonic_info))
                    changes.append((einfo, 0))
            einfo, score = sorted(changes, key=lambda l:l[1], reverse=True)[0]
            if score < 0.5 or not success:
                warnings.warn("Poor overlapping ({:.3f}) {}".format(score, l.strip()))
                continue
            einfo = "{},".format(einfo)
            fields[2] = einfo
            print('\t'.join(fields) + '\t{}'.format(score))

