#!/usr/bin/env python3

import argparse

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    with open(args.vcf) as infile:
        for l in infile:
            if l.startswith("#"):
                continue
            chrom, pos, name, ref, alt = l.strip().split('\t')[0:5]
            shift = 0
            while ref[shift] == alt[shift]:
                shift += 1
                if shift == len(ref) or shift == len(alt):
                    break
            ref, alt = ref[shift:], alt[shift:]
            pos = int(pos) + shift - 1
            print("{}\t{}\t{}\t{}".format(chrom, pos, pos + len(ref), name))

