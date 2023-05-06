#!/usr/bin/env python3

import argparse, os, sys, time
import json
from utils import copen

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument("--annovar", required=True, help="variant_function exonic_variant_function")
    p.add_argument('--start', required=True)
    #p.add_argument('--seed', type=int, default=2020)
    return p


def load_exonic_function(fn):
    exonic_function = dict()
    with copen(fn) as infile:
        for l in infile:
            nr, outcome, einfo, _, _, _, ref, alt = l.strip().split('\t')[0:8]
            nr = int(nr.replace("line", '')) - 1
            exonic_function[nr] = (outcome, einfo, ref, alt)
    return exonic_function


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)
    exonic_functions = load_exonic_function(args.annovar)
    start_change = {int(k.replace('line', '')) - 1: v for k, v in json.load(open(args.start)).items()}
    # metainfo = dict()
    print("##{}\n##{}\n#key\tlabel\tgroup\tchrom\tposition\tref\talt\tregion\tsplicing\tcoding\tcds_fs\tcds_indel\tstartgain\tstartloss\tstopgain\tstoploss\tsynonymous\tlength".format(time.asctime(), ' '.join(sys.argv)))
    with copen(args.vcf) as infile:
        cnt = 0
        for l in infile:
            if l.startswith('#'):
                continue
            chrom, pos, name, ref, alt = l.strip().split('\t')[:5]
            key, region = name.split("@@")[0], name.split("@@")[2]

            if cnt in start_change:
                if start_change[cnt] == "startloss":
                    startgain, startloss = 0, 1
                else:
                    startgain, startloss = 1, 0
            else:
                startgain, startloss = 0, 0

            if cnt in exonic_functions:
                outcome, einfo, ref2, alt2 = exonic_functions[cnt]
                # coding  cds_fs  cds_indel  startgain  startloss  stopgain  stoploss  synonymous
                if "frameshift" in outcome:
                    outcome = '\t'.join([str(x) for x in [1, 1, 0, startgain, startloss, 0, 0, 0]])
                elif "deletion" in outcome or "insertion" in outcome:
                    outcome = '\t'.join([str(x) for x in [1, 0, 1, startgain, startloss, 0, 0, 0]])
                elif "stopgain" in outcome:
                    outcome = '\t'.join([str(x) for x in [1, 0, 1, startgain, startloss, 1, 0, 0]])
                elif "stoploss" in outcome:
                    outcome = '\t'.join([str(x) for x in [1, 0, 1, startgain, startloss, 0, 1, 0]])
                elif outcome == "synonymous SNV":
                    outcome = '\t'.join([str(x) for x in [1, 0, 0, startgain, startloss, 0, 0, 1]])
                else:
                    outcome = '\t'.join([str(x) for x in [1, 0, 0, startgain, startloss, 0, 0, 0]])
            else:
                outcome = '\t'.join([str(0) for _ in range(8)])
            cnt += 1

            shift = 0
            while ref[shift] == alt[shift]:
                shift += 1
                if shift >= len(ref) or shift >= len(alt):
                    break
            ref_, alt_ = ref[shift:], alt[shift:]

            if len(ref_) == 1 and len(alt_) == 1:
                length = 0
            elif len(alt_) >= len(ref_):
                length = len(alt_)
            else:
                length = -len(ref_)
            if "splicing" in region:
                splicing = 1
            else:
                splicing = 0

            print('\t'.join([key, ".", ".", chrom, pos, ref, alt, region, str(splicing), outcome, str(length)]))
            # metainfo[key] = (label, group, function, outcome)


# frameshift deletion
# frameshift insertion
# frameshift substitution
# nonframeshift deletion
# nonframeshift insertion
# nonframeshift substitution
# nonsynonymous SNV
# # stopgain
# # stoploss
# # synonymous SNV
# unknown
