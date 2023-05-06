#!/usr/bin/env python3

import argparse
import os
import sys
from utils import copen

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('features', type=str)
    p.add_argument("-p", "--prefix", required=True)
    return p


if __name__ == "__main__":
    args = get_args().parse_args()

    header = list()

    ps = open(args.prefix + ".protein_snp.tsv", 'w')
    po = open(args.prefix + ".protein_other.tsv", 'w')
    ns = open(args.prefix + ".nonpro_snp.tsv", 'w')
    no = open(args.prefix + ".nonpro_other.tsv", 'w')
    ps.write("##protein SNP\n")
    po.write("##protein other\n")
    ns.write("##protein SNP\n")
    no.write("##protein other\n")

    with copen(args.features) as infile:
        for l in infile:
            if l.startswith("#"):
                ps.write(l)
                po.write(l)
                ns.write(l)
                no.write(l)
                if l.startswith("#key"):
                    fields = l.strip().split('\t')
                    length_col = fields.index("length")
                    prot_col = fields.index("PSSM-change.sumlog")
            else:
                fields = l.strip().split('\t')
                if fields[length_col] == '0' and fields[prot_col] != "NaN":
                    ps.write(l)
                elif fields[length_col] != '0' and fields[prot_col] != "NaN":
                    po.write(l)
                elif fields[length_col] == '0' and fields[prot_col] == "NaN":
                    ns.write(l)
                elif fields[length_col] != '0' and fields[prot_col] == "NaN":
                    no.write(l)

    ps.close()
    po.close()
    ns.close()
    no.close()

