#!/usr/bin/env python3

import argparse, os, sys, time, tqdm
import numpy as np
from utils import copen

import pyBigWig


def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('--bed', required=True)
    p.add_argument('-bw', nargs='+', required=True)
    return p

def fetch_signal(bed, bigwig, name):
    samples = list()
    with copen(bed) as infile:
        for l in infile:
            chrom, start, end, name = l.strip().split('\t')[0:4]
            key = name.split('@@')[0]
            samples.append((chrom, int(start), int(end), key))
    bw = pyBigWig.open(bigwig)
    sample_signal = dict()
    for chrom, start, end, key in tqdm.tqdm(samples, desc=os.path.basename(bigwig)):
        if chrom not in bw.chroms():
            signals1   = np.array([0], dtype=np.float32)
            signals10  = np.array([0], dtype=np.float32)
            signals100 = np.array([0], dtype=np.float32)
        else:
            try:
                signals1   = np.nan_to_num(np.array(bw.values(chrom, max(0, start),       min(bw.chroms()[chrom], max(start + 1, end))), dtype=np.float32))
                signals10  = np.nan_to_num(np.array(bw.values(chrom, max(0, start) - 10,  min(end + 10 ,          bw.chroms()[chrom])), dtype=np.float32))
                signals100 = np.nan_to_num(np.array(bw.values(chrom, max(0, start) - 100, min(end + 100,          bw.chroms()[chrom])), dtype=np.float32))
            except:
                raise RuntimeError("{}".format((bigwig, chrom, start, end)))
        sample_signal[key] = (signals1, signals10, signals100)
    return name, sample_signal


def signed_max(ar):
    if np.max(ar) == np.max(np.abs(ar)):
        m = np.max(ar)
    else:
        m = np.min(ar)
    return m

if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)

    sample_signals = dict()
    assay_list = list()
    for fn in sorted(args.bw):
        eid, assay = os.path.basename(fn).replace('.pval.signal.bigwig', '').split('-')
        sample_signals[assay] = fetch_signal(args.bed, fn, assay)[1]
        assay_list.append(assay)

    print("##{}\n##{}\n#key\t{}".format(time.asctime(), ' '.join(sys.argv), '\t'.join([
            "{}-mean\t{}-max\t{}.10bp-mean\t{}.10bp-max\t{}.100bp-mean\t{}.100bp-max".format(a, a, a, a, a, a) for a in assay_list
        ])), end='')
    with copen(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('@@')[0]
            print('\n' + key, end='')
            for k in assay_list:
                signals = sample_signals[k][key]
                print("\t{:.3f}\t{:.3f}".format(np.arcsinh(np.mean(signals[0])), np.arcsinh(signed_max(signals[0]))), end='')
                print("\t{:.3f}\t{:.3f}".format(np.arcsinh(np.mean(signals[1])), np.arcsinh(signed_max(signals[1]))), end='')
                print("\t{:.3f}\t{:.3f}".format(np.arcsinh(np.mean(signals[2])), np.arcsinh(signed_max(signals[2]))), end='')

            # print("{}\t{}\t{}".format(key, len(var_rbp_count["K562"][key]), len(var_rbp_count["HepG2"][key])))


