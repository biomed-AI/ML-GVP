#!/usr/bin/env python3

import argparse
import json
import os
import sys
import numpy as np
import pandas as pd
import gzip
import re

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("coding_change", help="marked coding change fasta file")
    # p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    args = get_args().parse_args()

    out = dict()
    re_line = re.compile(r"line(\d+)\s")
    with gzip.open(args.coding_change, 'rt') as infile:
        for l in infile:
            if l.startswith(">"):
                var_type = l.strip().split('|')[-1]
                if var_type == "startloss" or var_type == "startgain":
                    ans = re_line.findall(l)
                    line = "line" + ans[0]
                    out[line] = var_type
    
    json.dump(out, fp=sys.stdout, indent=4)
