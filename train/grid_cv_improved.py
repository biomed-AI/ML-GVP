#!/home/dingml/anaconda3/envs/xgboost/bin/python

import pickle
import argparse, os, sys, time, copy
import numpy as np
import json
import pandas as pd
from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score
from sklearn.model_selection import cross_val_score, GroupKFold, RandomizedSearchCV
import xgboost as xgb
from utils import copen, fetchFeaturesFromYaml
from skopt import BayesSearchCV
from skopt.space import Real, Integer

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
            data.append(pd.read_csv(fn, delimiter='\t', index_col=0, skiprows=n))
            pickle.dump(data[-1], open(cache, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    data = pd.concat(data, axis=0)
    label = data[args.label].to_numpy()
    groups = data[args.g].to_numpy()

    # feature_list = sorted(json.load(open(args.features)))
    feature_list = sorted(fetchFeaturesFromYaml(args.features))
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

    clf = BayesSearchCV(
            estimator=xgb_model,
            search_spaces={
                    "n_estimators": Integer(20, 300),
                    "max_depth": Integer(1, 10),
                    "gamma": [0, 0.001, 0.01, 0.1, 0.5, 1],
                    "min_child_weight": Integer(2, 10),
                    "subsample": Real(0.7, 1),
                    "learning_rate": [0.1, 0.3, 0.5, 0.7]
                },
            scoring=scoring,
            n_iter=100,
            random_state=args.seed,
            cv=splitter,
            verbose=5,
            n_jobs=3,
        )
    
    search = clf.fit(features, label)
    print("# best_score: {}".format(search.best_score_))
    print("# best_params: {}".format(search.best_params_))
    print("# best_results: {}".format(search.optimizer_results_))
    print("# scorer_: {}".format(search.scorer_))
    print("# model: {}".format(search.best_estimator_))
    # print("# params: ")
    json.dump(search.best_estimator_.get_params(), fp=open(args.o, 'w'), indent=4)

    # print(cross_val_score(estimator=xgb_model, X=features, y=labels, groups=groups, cv=splitter, scoring="average_precision"))
