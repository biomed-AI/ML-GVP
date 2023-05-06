#!/usr/bin/env python3

import argparse, os, sys, time

from collections import defaultdict
from utils import copen

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('--rbp', required=True, nargs='+')

    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)

    var_rbp_count = dict()
    for fn in args.rbp:
        cell, rbp = os.path.basename(fn).replace(".tsv.gz", '').split('_')[-2:]
        if cell not in var_rbp_count:
            var_rbp_count[cell] = defaultdict(set)
        with copen(fn) as infile:
            for l in infile:
                if not l.startswith("chr"):
                    continue
                name = l.split('\t')[3]
                key = name.split('@@')[0]
                var_rbp_count[cell][key].add(rbp)

    print("##{}\n##{}\n#key\tHepG2-RBP\tK562-RBP".format(time.asctime(), ' '.join(sys.argv)))
    with copen(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('@@')[0]
            print("{}\t{}\t{}".format(key, len(var_rbp_count["K562"][key]), len(var_rbp_count["HepG2"][key])))


