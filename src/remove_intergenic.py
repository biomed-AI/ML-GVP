#!/home/dingml/anaconda3/envs/xgboost/bin/python

import os, sys
from webbrowser import get
import pandas as pd
import json
import argparse
import numpy as np

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-r', '--ref', required=True)
    p.add_argument('-o', "--out", required=True)
    return p

def count_annotations(path):
    cnt = 0
    with open(path, mode='r') as infile:
        for l in infile:
            if l.startswith('##'):
                cnt += 1
            else:
                break
        return cnt


if __name__ == "__main__":
    args = get_args().parse_args()
    skip_rows = count_annotations(path=args.vcf)
    df = pd.read_csv(args.vcf, sep='\t', skiprows=skip_rows, header=0)
    np_ref = list()
    with open(args.ref, mode='r') as infile:
        for l in infile:
            np_ref.append(int(l.strip().rstrip('\n'))-1)  
    np_ref = np.array(np_ref, dtype=int)
    out_df = df.iloc[np_ref, :]
    out_df.to_csv(args.out, sep='\t', header=df.columns, index=None)