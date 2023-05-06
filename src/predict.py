#!/home/dingml/anaconda3/envs/xgboost/bin/python


import argparse, os, sys, time
import numpy as np
from typing import Any, Dict, List, Union
from grid_cv_improved import count_lines_begin_with

import logging, warnings, json, gzip, pickle
import pandas as pd
from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, f1_score
from sklearn.model_selection import cross_val_score, GroupKFold, RandomizedSearchCV
import xgboost as xgb
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
import re
# from biock.idle_gpu import idle_gpu

os.environ["CUDA_VISIBLE_DEVICES"] = '0' #str(idle_gpu())

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', "--input", nargs='+', required=True)
    p.add_argument('-m', required=True, help="model")
    p.add_argument('--gpu', default=-1, type=int)
    p.add_argument('-l', default="label")
    p.add_argument('-p', required=True)
    p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    np.random.seed(args.seed)

    if args.gpu >= 0:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)

    df = None
    for fn in args.input:
        if df is None:
            n_skip = count_lines_begin_with(fn, char="##")
            df = pd.read_csv(fn, sep='\t', index_col=0, skiprows=n_skip) 
        else:
            df = pd.concat((df, pd.read_csv(fn, sep='\t', index_col=0)), axis=0)
    print("##Input shape: {}".format(df.shape))

    model, feat_names = pickle.load(open(args.m, 'rb'))
    feat_importance = sorted(list(zip(feat_names, model.feature_importances_)), key=lambda x:x[1], reverse=True)
    d = dict()
    for f, s in feat_importance:
        d[f] = float(s)
    # json.dump(d, open("{}.feature_importances.json".format(args.m), 'w'), indent=4)

    keys = np.array(df.index, dtype=str).reshape(-1)
    features = df[feat_names].replace('.', 0).fillna(value=0).astype(np.float32)
    # features = (features - features.min(axis=0))/(features.max(axis=0) - features.min(axis=0))
    labels = np.array(df[args.l]).reshape(-1)

    print("## features: {}".format(features.shape))
    # print("## feature names: {}".format(feat_names))
    print("## labels: {}".format(np.unique(labels, return_counts=True)))

    prob = model.predict_proba(features).T[1]
    with open("./{}.prediction.txt".format(args.p), 'w') as out:
        for y, p, key in map(list, zip(labels, prob, keys)):
            out.write("{}\t{}\t{}\n".format(
                # y, '\t'.join(["{:.3f}".format(x) for x in prob]), key
                y, "{:.5f}".format(p), key
                ))

    if len(np.unique(labels)) > 1:
        try:
            labels = labels.astype(int)
            print("AUC/auPRC/F1:\t{:.4f}\t{:.4f}\t{:.4f}".format(roc_auc_score(labels, prob), \
                                                                 average_precision_score(labels, prob), \
                                                                 f1_score(labels, (prob> 0.5).astype(int))))
        except:
            exit(0)

    

