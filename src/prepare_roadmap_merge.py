#!/usr/bin/env python3

import argparse, os, sys, time
import numpy as np
from collections import defaultdict
from utils import copen


def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-a', '--roadmap-annotations', dest="anno", nargs='+', required=True)
    #p.add_argument('--seed', type=int, default=2020)
    return p


def load_scores(fn):
    scores = dict()
    with copen(fn) as infile:
        for l in infile:
            if l.startswith('##'):
                continue
            elif l.startswith('#key'):
                names = l.strip().split('\t')[1:]
            else:
                fields = l.strip().split('\t')
                key = fields[0]
                scores[key] = dict()
                for c, v in map(list, zip(names, fields[1:])):
                    scores[key][c] = float(v)
    return scores


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    assay_scores = dict()
    assay_list = set()
    for fn in args.anno:
        sample_signals = load_scores(fn)
        for k in sample_signals:
            if k not in assay_scores:
                assay_scores[k] = defaultdict(list)
            for c, v in sample_signals[k].items():
                assay_list.add(c)
                assay_scores[k][c].append(v)
    # json.dump(assay_scores, fp=open("out.json", 'w'), indent=4)

    assay_list = sorted(list(assay_list))
    print("##{}\n##{}\n#key\t{}".format(time.asctime(), ' '.join(sys.argv), '\t'.join(assay_list)), end='')
    with copen(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('@@')[0]
            print('\n' + key, end='')
            for a in assay_list:
                if a in assay_scores[key]:
                    signals = assay_scores[key][a]
                    print("\t{:.3f}".format(np.max(signals)), end='')
                else:
                    print("\t0", end='')

            # print("{}\t{}\t{}".format(key, len(var_rbp_count["K562"][key]), len(var_rbp_count["HepG2"][key])))


