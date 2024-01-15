#!/usr/bin/env python3

import argparse, os, sys, time
import pandas as pd
import gzip
import json
from tqdm import tqdm
import time
import numpy as np
import xgboost as xgb
import warnings
warnings.filterwarnings('ignore')
import tabix
from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, f1_score

LABELS = {
    'Benign': 0,
    'Benign/Likely_benign': 0,
    'Likely_pathogenic': 1,
    'risk_factor': 1,
    'Pathogenic/Likely_pathogenic': 1,
    'drug_response': 1,
    'Uncertain_significance': 0,
    'Likely_benign': 0,
    'Pathogenic': 1,
    'Conflicting_interpretations_of_pathogenicity': 0,
}

# LABELS = {"LOF": 1,
#             "FUNC": 0,
#             "INT": 0}

# LABELS = {
#     "wt-like": 1,
#     "possibly_wt-like":1,
#     "possibly_low": 0,
#     "low": 0}



def getargs():
    args = argparse.ArgumentParser()
    args.add_argument("--ref_dir", required=True, help='the reference BED file for all predictors(directory)')
    args.add_argument("-i", required=True, help="the input variants file")
    return args.parse_args()

if __name__ == "__main__":
    args = getargs()

    # CADD for variants
    method = "CADD"
    print(f"{'='*30}{method}-BEGIN{'='*30}")
    tb_snp = tabix.open(f"{args.ref_dir}/whole_genome_SNVs_CADD-PHRED.tsv.gz")
    tb_indel = tabix.open(f"{args.ref_dir}/InDels_inclAnno.tsv.gz")
    fn = args.i
    cadd_score = dict()
    with open(fn, mode='r') as infile:
        for line in tqdm(infile):
            if line.startswith("#"):
                continue
            fields = line.strip("\b").split("\t")
            key = fields[0]
            chrom, pos, ref, alt = key.split("_")
            score = -1.0
            if len(ref) == 1 and len(alt) == 1:
                records = tb_snp.query(chrom[3:], int(pos)-1, int(pos)) # 0-based
                for record in records:
                    if record[2] == ref and record[3] == alt:
                        score = float(record[-1])
            else:
                records = tb_indel.query(chrom[3:], int(pos)-1, int(pos)) # 0-based
                for record in records:
                    if record[2] == ref and record[3] == alt:
                        score = float(record[-1])
            cadd_score[key] = score
    
    print(f"{method} scores were save at {args.i}.{method}")
    with open(f"{args.i}.{method}", mode='w') as outfile:
        for _, key in enumerate(cadd_score.keys()):
            outfile.writelines("{}\t{}\n".format(key, cadd_score[key]))
    print(f"{'='*30}{method}-END{'='*30}")


    # MVP for missense
    method = "MVP"
    print(f"{'='*30}{method}-BEGIN{'='*30}")
    tb = tabix.open(f"{args.ref_dir}/MVP_score_hg19.vcf.gz")
    fn = args.i
    mvp_score = dict()
    with open(fn, mode='r') as infile:
        for line in tqdm(infile):
            if line.startswith("#"):
                continue
            fields = line.strip("\b").split("\t")
            key = fields[0]
            chrom, pos, ref, alt = key.split("_")
            records = tb.query(chrom, int(pos)-1, int(pos)) # 0-based
            score = -1.0
            for record in records:
                if record[3] == ref and record[4] == alt:
                    score = float(record[-1].split("|")[-1].split(":")[-1])
            mvp_score[key] = score
    
    print(f"{method} scores were save at {args.i}.{method}")
    with open(f"{args.i}.{method}", mode='w') as outfile:
        for _, key in enumerate(mvp_score.keys()):
            outfile.writelines("{}\t{}\n".format(key, mvp_score[key]))
    print(f"{'='*30}{method}-END{'='*30}")

    # M-CAP for missense
    method = "MCAP"
    print(f"{'='*30}{method}-BEGIN{'='*30}")
    tb = tabix.open(f"{args.ref_dir}/mcap_v1_4.vcf.gz")
    fn = args.i
    mcap_score = dict()
    with open(fn, mode='r') as infile:
        for line in tqdm(infile):
            if line.startswith("#"):
                continue
            fields = line.strip("\b").split("\t")
            key = fields[0]
            chrom, pos, ref, alt = key.split("_")
            records = tb.query(chrom, int(pos)-1, int(pos)) # 0-based
            score = -1.0
            for record in records:
                if record[3] == ref and record[4] == alt:
                    score = float(record[-1].split("|")[0].split(":")[-1])
            mcap_score[key] = score
    
    print(f"{method} scores were save at {args.i}.{method}")
    with open(f"{args.i}.{method}", mode='w') as outfile:
        for _, key in enumerate(mcap_score.keys()):
            outfile.writelines("{}\t{}\n".format(key, mcap_score[key]))
    print(f"{'='*30}{method}-END{'='*30}")


    # PrimateAI for missense
    method = "PrimateAI"
    print(f"{'='*30}{method}-BEGIN{'='*30}")
    tb = tabix.open(f"{args.ref_dir}/PrimateAI_scores_v0.2.vcf.gz")
    fn = args.i
    primateai_score = dict()
    with open(fn, mode='r') as infile:
        for line in tqdm(infile):
            if line.startswith("#"):
                continue
            fields = line.strip("\b").split("\t")
            key = fields[0]
            chrom, pos, ref, alt = key.split("_")
            records = tb.query(chrom, int(pos)-1, int(pos)) # 0-based
            score = -1.0
            for record in records:
                if record[3] == ref and record[4] == alt:
                    score = float(record[-1].split("|")[-1].split(":")[-1])
            primateai_score[key] = score
    
    print(f"{method} scores were save at {args.i}.{method}")
    with open(f"{args.i}.{method}", mode='w') as outfile:
        for _, key in enumerate(primateai_score.keys()):
            outfile.writelines("{}\t{}\n".format(key, primateai_score[key]))
    print(f"{'='*30}{method}-END{'='*30}")

    # AlphaMissense for missense
    method = "AlphaMissense"
    print(f"{'='*30}{method}-BEGIN{'='*30}")
    tb = tabix.open(f"{args.ref_dir}/AlphaMissense_hg19.vcf.gz")
    fn = args.i
    alphamissense_score = dict()
    with open(fn, mode='r') as infile:
        for line in tqdm(infile):
            if line.startswith("#"):
                continue
            fields = line.strip("\b").split("\t")
            key = fields[0]
            chrom, pos, ref, alt = key.split("_")
            records = tb.query(chrom, int(pos)-1, int(pos)) # 0-based
            score = -1.0
            for record in records:
                if record[3] == ref and record[4] == alt:
                    score = float(record[-1].split("|")[-2].split(":")[-1])
            alphamissense_score[key] = score
    
    print(f"{method} scores were save at {args.i}.{method}")
    with open(f"{args.i}.{method}", mode='w') as outfile:
        for _, key in enumerate(alphamissense_score.keys()):
            outfile.writelines("{}\t{}\n".format(key, alphamissense_score[key]))
    print(f"{'='*30}{method}-END{'='*30}")







