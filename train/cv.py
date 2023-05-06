#!/home/dingml/anaconda3/envs/xgboost/bin/python

import argparse, os, sys, time
from textwrap import indent
import numpy as np
from typing import Any, Dict, List, Union, final
from functools import partial
print = partial(print, flush=True)

import json, pickle
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, f1_score
import xgboost as xgb
from grid_cv_improved import count_lines_begin_with
from utils import fetchFeaturesFromYaml


def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-d', "--data", nargs='+', required=True)
    p.add_argument('-c', required=True, help="config")
    p.add_argument('-l', "--label", default="label")
    p.add_argument('-g', default="chrom", help="group")
    p.add_argument('-f', "--features", required=True, help="feature list (yml)")
    p.add_argument("--ignored", help="igored features", nargs='+')
    p.add_argument('-p', required=True, help="prefix")
    p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    np.random.seed(args.seed)

    print("##{}\n##{}".format(time.asctime(), ' '.join(sys.argv)))
    data = list()
    for fn in args.data:
        print(f"## loading {fn}")
        cache = fn + ".cache.pkl"
        if os.path.exists(cache):
            data.append(pickle.load(open(cache, 'rb')))
        else:
            n = count_lines_begin_with(fn, char="##")
            data.append(pd.read_csv(fn, sep='\t', index_col=0, skiprows=n))
            pickle.dump(data[-1], open(cache, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    data = pd.concat(data, axis=0)
    label = data[args.label].to_numpy()
    groups = data[args.g].to_numpy()

    # feature_list = sorted(json.load(open(args.features)))
    
    feature_list = sorted(fetchFeaturesFromYaml(args.features, ignored_features=[] if args.ignored==None else args.ignored))
    features = data[feature_list].fillna(value=0).replace('.', 0).astype(np.float32)
    # features = (features - features.min(axis=0))/(features.max(axis=0) - features.min(axis=0))

    print("## features: {}".format(features.shape))
    print("## feature names: {}".format(feature_list))
    print("## groups: {}".format(np.unique(groups, return_counts=True)))
    print("## labels: {}".format(np.unique(label, return_counts=True)))
    keys = features.index.to_numpy()

    splitter = GroupKFold(n_splits=5)

    config = json.load(open(args.c))
    config["n_jobs"] = 4

    final_label, final_prob, final_sample = list(), list(), list()
    for fold, (train_idx, valid_idx) in enumerate(splitter.split(features, groups=groups)):
        train_X = features.iloc[train_idx, :]
        train_y = label[train_idx]
        xgb_model = xgb.XGBClassifier(**config)
        xgb_model.fit(train_X, train_y)
        prob = xgb_model.predict_proba(features.iloc[valid_idx, :])
        final_label.append(label[valid_idx])
        final_prob.append(prob)
        final_sample.append(keys[valid_idx])

    final_label = np.concatenate(final_label)
    final_prob = np.concatenate(final_prob, axis=0).T[1]
    final_sample = np.concatenate(final_sample)
    
    with open("./predictions/{}.prediction.txt".format(args.p), 'w') as out:
        for y, prob, key in map(list, zip(final_label, final_prob, final_sample)):
            out.write("{}\t{}\t{}\n".format(
                # y, '\t'.join(["{:.3f}".format(x) for x in prob]), key
                y, "{:.5f}".format(prob), key
                ))
    # if len(np.unique(label)) > 1:
    #     try:
    #         labels = label.astype(int)
    #         print("AUC/auPRC/F1:\t{:.4f}\t{:.4f}\t{:.4f}".format(roc_auc_score(final_label, final_prob), \
    #                                                              average_precision_score(final_label, final_prob), \
    #                                                              f1_score(final_label, (final_prob> 0.5).astype(int))))
    #     except:
    #         exit(0)

    final_model = xgb.XGBClassifier(**config)
    final_model.fit(features, label)

    pickle.dump((final_model, feature_list), open("../models/{}.model.pkl".format(args.p), 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

