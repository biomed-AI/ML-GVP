#!/usr/bin/env python3

import argparse
import os
import sys
from utils import copen

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("prediction")
    # p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    args = get_args().parse_args()

    cache = list()
    with open(args.prediction) as infile:
        for l in infile:
            label, prob, key = l.strip().split('\t')
            idx, key = key.split("@")
            cache.append((int(idx), key, label, prob))
    
    cache = sorted(cache, key=lambda x:x[0])
    print("#key\tlabel\tscore")
    for _, key, label, prob in cache:
        print(f"{key}\t{label}\t{prob}")
