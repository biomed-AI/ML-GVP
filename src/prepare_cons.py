#!/usr/bin/env python3

import argparse, os, sys, time, tqdm
import numpy as np

from collections import OrderedDict
import pyBigWig
from utils import copen


def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('--bed', required=True)
    p.add_argument("--name", required=True, nargs='+')
    p.add_argument("--bw", required=True, nargs='+')

    #p.add_argument('--seed', type=int, default=2020)
    return p

def fetch_signal(bed, bigwig, name):
    samples = list()
    with copen(args.bed) as infile:
        for l in infile:
            chrom, start, end, name = l.strip().split('\t')[0:4]
            key = name.split('@@')[0]
            samples.append((chrom, int(start), int(end), key))
    bw = pyBigWig.open(bigwig)
    sample_signal = dict()
    for chrom, start, end, key in tqdm.tqdm(samples, desc=os.path.basename(bigwig)):
        if chrom not in bw.chroms():
            signal1 = np.array([0], dtype=np.float32)
            signal10 = np.array([0], dtype=np.float32)
            signal100 = np.array([0], dtype=np.float32)
        else:
            try:
                signal1   = np.nan_to_num(np.array(bw.values(chrom, start, max(start + 1, end)), dtype=np.float32))
            except:
                raise RuntimeError("{}".format((bigwig, chrom, start, end)))
            try:
                signal10  = np.nan_to_num(np.array(bw.values(chrom, start - 10,  end + 10), dtype=np.float32))
            except:
                raise RuntimeError("{}".format((bigwig, chrom, start, end, 10)))
            try:
                signal100 = np.nan_to_num(np.array(bw.values(chrom, start - 100, end + 100), dtype=np.float32))
            except:
                raise RuntimeError("{}".format((bigwig, chrom, start, end, 100)))
        sample_signal[key] = (signal1, signal10, signal100)
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
    
    # job_list = [(args.bed, v, k) for (k, v) in cons_files.items()]
    # with Pool(processes=2) as pool:
    #     res = pool.starmap(fetch_signal, job_list)
    # cons_files = OrderedDict({
    #     "phastCons100way": "/home/chenken/db/UCSC/hg19/conservation/hg19.100way.phastCons.bw",
    #     "phyloP100way": "/home/chenken/db/UCSC/hg19/conservation/hg19.100way.phyloP100way.bw",
    #     "phastCons46way.primates": "/home/chenken/db/UCSC/hg19/conservation/primates.phastCons46way.bw",
    #     "phyloP46way.primates": "/home/chenken/db/UCSC/hg19/conservation/primates.phyloP46way.bw",
    #     "phastCons46way.vertebrate": "/home/chenken/db/UCSC/hg19/conservation/vertebrate.phastCons46way.bw",
    #     "phyloP46way.vertebrate": "/home/chenken/db/UCSC/hg19/conservation/vertebrate.phyloP46way.bw"
    #     })
    assert len(args.name) == len(args.bw), "{}".format(((len(args.name), args.name), (len(args.bw), args.bw)))
    cons_files = OrderedDict()
    for name, bw in zip(args.name, args.bw):
        cons_files[name] = bw
    sample_signals = dict()
    for k, v in cons_files.items():
        sample_signals[k] = fetch_signal(args.bed, v, k)[1]

    print("##{}\n##{}\n#key\t{}".format(time.asctime(), ' '.join(sys.argv), '\t'.join([
            "{}-mean\t{}-max\t{}.10bp-mean\t{}.10bp-max\t{}.100bp-mean\t{}.100bp-max".format(k, k, k, k, k, k) for k in cons_files
        ])), end='')
    with copen(args.vcf) as infile:
        for l in infile:
            if l.startswith('#'):
                continue
            key = l.split('\t')[2].split('@@')[0]
            print('\n' + key, end='')
            for k in cons_files:
                signals = sample_signals[k][key]
                print("\t{:.3f}\t{:.3f}".format(np.mean(signals[0]), signed_max(signals[0])), end='')
                print("\t{:.3f}\t{:.3f}".format(np.mean(signals[1]), signed_max(signals[1])), end='')
                print("\t{:.3f}\t{:.3f}".format(np.mean(signals[2]), signed_max(signals[2])), end='')

            # print("{}\t{}\t{}".format(key, len(var_rbp_count["K562"][key]), len(var_rbp_count["HepG2"][key])))


