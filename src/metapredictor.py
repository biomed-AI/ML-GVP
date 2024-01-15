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
import pickle
import torch.nn.functional as F


def getargs():
    args = argparse.ArgumentParser()
    args.add_argument("-i", required=True, help="The input variants file")
    args.add_argument("-o", required=True, help="The output txt file")
    args.add_argument("-m", required=True, help="The meta model path")
    return args.parse_args()


if __name__ == "__main__":
    # Test with ensemble model
    args = getargs()
    fn = args.i
    savepath = args.m
    print(f"## loading {fn}")
    print(f"## loading model {savepath}")
    model, feature_list = pickle.load(open(savepath, 'rb'))
    data = pd.read_csv(fn, delimiter="\t", index_col=0)
    keys = data.index
    label  = np.array(data['label'])
    features = np.array(data[feature_list])
    print("## labels: {}".format(np.unique(label, return_counts=True)))
    print("## feature name: {}-{}".format(feature_list, features.shape))

    prob = model.predict_proba(features).T[1]
    prob = prob/np.max(prob)
    if len(np.unique(label)) != 1:
        print(
            roc_auc_score(y_true=label, y_score=prob),
            average_precision_score(y_true=label, y_score=prob)
        )
    ensemble_dict = dict()
    with open(args.o, mode='w') as outfile:
        outfile.writelines("#key\tmeta-score\n")
        for i, key in enumerate(keys):
            outfile.writelines(f"{key}\t{str(prob[i])}\n")
