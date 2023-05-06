#!/usr/bin/env python3

import argparse, os, sys, time

from collections import defaultdict
from utils import copen

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-isec', nargs='+', required=True)

    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    sample_peaks = dict()
    assay_list = set()
    for fn in args.isec:
        bn = os.path.basename(fn)
        assay = bn.replace(".tsv.gz", '').replace(".tsv", "").split('_')[-1]
        assay_list.add(assay)
        with copen(fn) as infile:
            for l in infile:
                if l.startswith('#'):
                    continue
                key = l.split('\t')[3].split('@@')[0]
                if key not in sample_peaks:
                    sample_peaks[key] = defaultdict(set)
                sample_peaks[key][assay].add(bn)

    assay_list = sorted(list(assay_list))
    print("##{}\n##{}\n#key\t{}".format(time.asctime(), ' '.join(sys.argv), '\t'.join(["{}".format(a) for a in assay_list])))
    with open(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('@@')[0]
            if key in sample_peaks:
                print('{}\t{}'.format(key, '\t'.join([str(len(sample_peaks[key][a])) for a in assay_list])))
            else:
                print('{}\t{}'.format(key, '\t'.join([str(0) for a in assay_list])))
            # print("{}\t{}\t{}".format(key, len(var_rbp_count["K562"][key]), len(var_rbp_count["HepG2"][key])))


