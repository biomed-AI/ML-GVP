#!/usr/bin/env python3

import argparse, re
import warnings, gzip
from collections import OrderedDict


patterns = OrderedDict({
    "point": re.compile(r'protein-altering  \(position (\d+) changed from ([A-Z]) to ([A-Z\*]+)?\)'),
    "indel": re.compile(r'protein-altering  \(position (\d+)-(\d+) changed from ([A-Z\*]*) to ([A-Z\*]*)\)'),
    "stopgain-point": re.compile(r'immediate-stopgain  \(position (\d+) changed from ([A-Z\*]*) to ([A-Z]*)\)'),
    "stopgain-range": re.compile(r'immediate-stopgain  \(position (\d+)-(\d+) changed from ([A-Z\*]*) to ([A-Z\*]*)\)'),
    "stoploss-point": re.compile(r'immediate-stoploss  \(position (\d+) changed from ([A-Z\*]+) to ([A-Z\*]*)\)'),
    "stoploss-range": re.compile(r'immediate-stoploss  \(position (\d+)-(\d+) changed from ([A-Z\*]+) to ([A-Z\*]*)\)'),
    "insertion": re.compile(r"protein-altering  \(position (\d+)-(\d+) has insertion ([A-Z\*]+)\)"),
    "startloss": re.compile(r'startloss'),
    # "silent": re.compile(r'p.([A-Z\*])(\d+)([A-Z\*]) silent  \(no amino acid change\)'),
    "silent": re.compile(r'silent  \(no amino acid change\)'),
})


def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('coding_change')
    p.add_argument('-vcf', required=True)
    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)


    raw_vcf = list()
    with open(args.vcf) as infile:
        for l in infile:
            if l.startswith("#"):
                continue
            raw_vcf.append(l.strip().split('\t')[0:5])

    name = None
    seq = list()
    # ref_seqs = dict()
    fasta = dict()
    altered_protein = list()
    with gzip.open(args.coding_change, 'rt') as infile:
        for l in infile:
            if l.startswith('>'):
                if name is not None:
                    if name.endswith("WILDTYPE"):
                        tx_id = name.split(' ')[1]
                        fasta[tx_id] = ''.join(seq)
                        # ref_seqs[tx_id] = ''.join(seq)
                    else:
                        fasta[name] = ''.join(seq)
                        altered_protein.append(name)
                    seq = list()
                name = l.strip().lstrip('>')
            else:
                seq.append(l.strip())

    if name.endswith("WILDTYPE"):
        tx_id = name.split(' ')[1]
        fasta[tx_id] = ''.join(seq)
        # ref_seqs[tx_id] = ''.join(seq)
    else:
        fasta[name] = ''.join(seq)
        altered_protein.append(name)
    
    for name in altered_protein:
        tx_id = name.split(' ')[1]
        raw_seq = fasta[tx_id]
        new_seq = fasta[name]
        found = False
        for vartype, pat in patterns.items():
            res = pat.findall(name)
            if len(res) == 1:
                if vartype == "point":
                    pos, ref, alt = res[0]
                    pos = int(pos)
                    left_len, right_len = pos - 1, len(raw_seq) - pos
                    if len(ref) == len(alt) and len(ref) == 1:
                        group = "missense_single"
                    else:
                        group = "indel"
                elif vartype == "stopgain-point" or vartype == "stoploss-point":
                    pos, ref, alt = res[0]
                    left_len = int(pos) - 1
                    right_len = len(raw_seq) - int(pos)
                    if vartype == "stopgain-point":
                        group = "stopgain_single"
                    else:
                        group = "stoploss_single"
                elif vartype == "indel" or vartype == "stopgain-range" or vartype == "stoploss-range":
                    start, end, ref, alt = res[0]
                    left_len, right_len = int(start) - 1, len(raw_seq) - int(end)
                    assert left_len + right_len + len(ref) == len(raw_seq)
                    if vartype == "indel":
                        group = "indel"
                    elif vartype == "stopgain-range":
                        group = "stopgain_range"
                    elif vartype == "stoploss-range":
                        group = "stoploss_range"
                elif vartype == "insertion":
                    start, end, ins = res[0]
                    left_len = int(start)
                    right_len = len(raw_seq) - left_len
                    ref, alt = "", ins
                    group = vartype
                elif vartype == "startloss":
                    left_len, right_len = 0, 0
                    ref = raw_seq
                    alt = new_seq
                    group = vartype
                elif vartype == "silent":
                    continue
                    # ref, pos, alt = res[0]
                    # pos = int(pos)
                    # left_len, right_len = pos - 1, len(raw_seq) - pos
                    # assert new_seq == raw_seq[:left_len] + alt + raw_seq[-right_len:]
                # if len(ref) > len(alt):
                #     alt_space = ' ' * (len(ref) - len(alt))
                #     ref_space = ""
                # elif len(alt) > len(ref):
                #     alt_space = ""
                #     ref_space = ' ' * (len(alt) - len(ref))
                # else:
                #     ref_space, alt_space = "", ""
                assert raw_seq[0:left_len] + alt + raw_seq[len(raw_seq)-right_len:] == new_seq, "{}\n{}\n{}\n{}".format((name, vartype, res), raw_seq, new_seq, raw_seq[0:left_len] + alt + raw_seq[len(raw_seq)-right_len:])
                nline = int(name.split()[0].replace("line", '')) - 1

                print(">{}|{}|{}\n{}[{}/{}]{}".format(
                    raw_vcf[nline][2].split('|')[0], name, group,
                    raw_seq[0:left_len], ref, alt, raw_seq[len(raw_seq)-right_len:]
                ))
                del group
                found = True
                break
        if not found:
            warnings.warn("{}".format(name))
    # json.dump(ref_seqs, open(args.r, 'w'), indent=4)
