#!/usr/bin/env python3

import argparse, os, sys, time
import numpy as np
from utils import copen

# import matplotlib.pyplot as plt
# from scipy.stats import pearsonr, spearmanr, ttest_ind, mannwhitneyu
# from sklearn.metrics import auc, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score

def pandas_df2dict(fn, delimiter='\t', **kwargs):
    import pandas as pd
    kwargs["delimiter"] = delimiter
    df = pd.read_csv(fn, **kwargs)
    d = dict()
    for k in df.columns:
        d[k] = np.array(df[k])
    return d

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-m', required=True, help="MMSPlice out")

    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    sample_mmsplice = dict()
    header = ["MMS_delta_logit_psi", "MMS_ref_acceptorIntron", "MMS_alt_acceptorIntron", "MMS_ref_acceptor", "MMS_alt_acceptor", "MMS_ref_exon", "MMS_alt_exon", "MMS_ref_donor", "MMS_alt_donor", "MMS_ref_donorIntron", "MMS_alt_donorIntron", "MMS_pathogenicity", "MMS_efficiency"]
    data = pandas_df2dict(args.m, delimiter=',', index_col=0)
    for i, mid in enumerate(data["ID"]):
        ar = [data[k.replace('MMS_', '')][i] for k in header]
        chrom, pos, ref_alt = mid.split(':')
        ref, alt = ref_alt.split('>')
        key = "{}_{}_{}_{}".format(chrom, pos, ref, alt)
        sample_mmsplice[key] = ar

    print("##{}\n##{}\n#key\t{}".format(time.asctime(), ' '.join(sys.argv), '\t'.join(["{}".format(a) for a in header])))
    with copen(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('@@')[0]
            mm_key = key.split('@')[1]
            if mm_key in sample_mmsplice:
                print('{}\t{}'.format(key, '\t'.join(["{:.3f}".format(x) for x in sample_mmsplice[mm_key]])))
            else:
                print('{}\t{}'.format(key, '\t'.join([str(0) for _ in header])))
            # print("{}\t{}\t{}".format(key, len(var_rbp_count["K562"][key]), len(var_rbp_count["HepG2"][key])))


