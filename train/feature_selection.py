#!/home/dingml/anaconda3/envs/xgboost/bin/python

import pickle
import argparse, os, sys, time ,copy
import numpy as np
import json
import pandas as pd
from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LassoCV
from sklearn.model_selection import GroupKFold
import xgboost as xgb
from utils import copen

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-d', "--data", nargs='+', required=True)
    p.add_argument("-l", "--label", default="label")
    p.add_argument('-o', required=True)
    p.add_argument('--device', default="cpu", type=str, help="cpu or gpu index")
    p.add_argument('-g', default="chrom", help="group key")
    p.add_argument('-f', "--features", default="./noncoding_features.json")
    p.add_argument('--seed', type=int, default=2020)
    return p

def count_lines_begin_with(fn, char, scan_all=False):
    cnt = 0
    with copen(fn) as infile:
        for l in infile:
            if l.startswith(char):
                cnt += 1
            elif not scan_all:
                break
    return cnt

if __name__ == "__main__":
    args = get_args().parse_args()
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
            data.append(pd.read_csv(fn, delimiter='\t', index_col=0, skiprows=n))
            pickle.dump(data[-1], open(cache, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    data = pd.concat(data, axis=0)
    label = data[args.label].to_numpy()
    groups = data[args.g].to_numpy()

    feature_list = sorted(json.load(open(args.features)))
    features = data[feature_list].fillna(value=0).replace('.', 0).astype(np.float32)

    print("## features: {}".format(features.shape))
    print("## groups: {}".format(np.unique(groups, return_counts=True)))
    print("## labels: {}".format(np.unique(label, return_counts=True)))

    splitter = GroupKFold(n_splits=5).split(X=features, groups=groups)
    config = dict(use_label_encoder=False, eval_metric="logloss", n_jobs=4)
    if args.device != "cpu":
        config["tree_method"] = 'gpu_hist'
    xgb_model = xgb.XGBClassifier(**config)

    scoring = "roc_auc"

    rfe = RFECV(
        estimator=xgb_model,
        step=1,
        min_features_to_select=1,
        cv=splitter,
        scoring=scoring,
        verbose=5,
        n_jobs=5,
    )

    selector = rfe.fit(features, label)
    
    
    print("feature support:{}".format(selector.support_))
    print("feature ranking:{}".format(selector.ranking_))
    feature_list = np.array(feature_list)
    with open(args.o, mode='w') as outfile:
        s = json.dumps(list(feature_list[np.where(selector.ranking_ == 1.0)[0]]), indent=4)
        outfile.write(s)