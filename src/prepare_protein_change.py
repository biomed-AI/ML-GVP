#!/home/dingml/anaconda3/envs/xgboost/bin/python

import argparse, os, sys, time
import logging, warnings, json, gzip, pickle
from utils import copen
# warning = logging.warning
from multiprocessing import Pool
from Bio.Align import PairwiseAligner, substitution_matrices

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
    p.add_argument('vcf')
    p.add_argument('-m', help="./final_train_092921.protein.mark_mut.tsv", required=True)
    
    #p.add_argument('--seed', type=int, default=2020)
    return p

def mp_align(name, ref_seq, alt_seq):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.mode = "global"
    ref_score = aligner.score(ref_seq, ref_seq)
    alt_score = aligner.score(ref_seq, alt_seq)
    score = (alt_score - ref_score) / ref_score
    return score
    # print("{}\t{}\t{}\t{:.5f}\t{:.5f}".format(name, ref, alt, score, 1 - score))
 


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    job_list = list()
    with gzip.open(args.m, 'rt') as infile:
        for l in infile:
            if l.startswith('>'):
                name = l.split('@@')[0].lstrip('>')
            else:
                left_ref, alt_right = l.strip().rstrip('*').split('/')
                left, ref = left_ref.split('[')
                alt, right = alt_right.split(']')
                left, ref, alt, right =left.rstrip("*"), ref.rstrip("*"), alt.rstrip("*"), right.rstrip("*")
                ref_seq = left + ref + right
                alt_seq = left + alt + right
                assert '*' not in ref_seq and '*' not in alt_seq, "{}".format(name)
                job_list.append((name, ref_seq, alt_seq))

                # score = aligner.score(ref_seq, alt_seq) / max(len(ref_seq), len(alt_seq)) / 2
                # print("{}\t{}\t{}\t{:.5f}\t{:.5f}".format(name, ref, alt, score, 1 - score))
    batch_size = 256

    n_batch = np.ceil(len(job_list) / batch_size).astype(int)



    var_prot_sim = dict()
    for i in range(n_batch):
        with Pool(processes=16) as pool:
            scores = pool.starmap(mp_align, job_list[i * batch_size : (i + 1) * batch_size])
    
        for j, score in enumerate(scores):
            name = job_list[i * batch_size + j][0]
            var_prot_sim[name] = score
    
    
    print("##{}\n##{}\n#key\tprot_sim_change".format(time.asctime(), ' '.join(sys.argv)))
    with open(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('@@')[0]
            if key in var_prot_sim:
                print('{}\t{:.5g}'.format(key, var_prot_sim[key]))
            else:
                print('{}\t{}'.format(key, 0))
            # print("{}\t{}\t{}".format(key, len(var_rbp_count["K562"][key]), len(var_rbp_count["HepG2"][key])))


