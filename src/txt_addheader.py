#!/usr/bin/env python3

import argparse, os, sys, time
import pandas as pd
import gzip
import json
from tqdm import tqdm
import time
import numpy as np
import xgboost as xgb
from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, f1_score



def getargs():
    args = argparse.ArgumentParser()
    args.add_argument("-i", required=True, help="The input txt file")
    args.add_argument("-o", required=True, help="The output txt file")
    return args.parse_args()


if __name__ == "__main__":
    args = getargs()
    fn = args.i
    with open(args.i, mode='r') as infile, open(args.o, mode='w') as outfile:
        l = f"#key\tlabel\tchrom\tCADD-score\tMVP-score\tM-CAP-score\tPrimateAI-score\tAlphaMissense-score\tML-GVE-score\n"
        outfile.writelines(l)
        for _, line in tqdm(enumerate(infile)):
            if line.startswith("#"):
                pass
            else:
                fields = line.strip("\n").split("\t")
                chrom = fields[0].split("_")[0][3:]
                if chrom == "X":
                    chrom = 23
                elif chrom == "Y":
                    chrom = 24
                else:
                    chrom = int(chrom)
                l = f"{fields[0]}\t{fields[1]}\t{chrom}\t{fields[2]}\t{fields[3]}\t{fields[4]}\t{fields[5]}\t{fields[6]}\t{fields[7]}\n"
                outfile.writelines(l)
