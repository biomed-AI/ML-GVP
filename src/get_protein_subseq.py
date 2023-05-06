#!/usr/bin/env python3

import argparse, gzip
#import pandas as pd
#import matplotlib.pyplot as plt
#from scipy.stats import pearsonr, spearmanr

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('marked_fa')
    p.add_argument("-pad", type=int, default=101)
    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    
    with gzip.open(args.marked_fa, 'rt') as infile:
        for l in infile:
            if l.startswith(">"):
                name = l.strip().lstrip('>')
            else:
                left_ref, alt_right = l.strip().split('/')
                left, ref = left_ref.split('[')
                alt, right = alt_right.split(']')
                left, ref, alt, right = left.strip('*'), ref.strip('*'), alt.strip('*'), right.strip('*')
                left_pad = min(max(20, (args.pad - len(ref)) // 2), len(left))
                right_pad = min(args.pad - len(ref) - left_pad, len(right))
                left_pad = min(args.pad - len(ref) - right_pad, len(left))
                if left_pad == 0 and right_pad == 0:
                    left_pad = min(len(left), 20)
                print(">{}|left/right_pad={}/{}|ref/alt={}/{}\n{}{}{}".format(name, left_pad, right_pad, ref, alt, left[-left_pad:], ref, right[:right_pad]))

