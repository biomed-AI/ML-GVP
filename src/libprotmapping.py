#!/usr/bin/env python3

import argparse, os, sys, re, tqdm
import warnings, gzip, pickle
# warning = logging.warning

import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#from scipy.stats import pearsonr, spearmanr
from typing import Any, Callable, Dict, List, OrderedDict, Union
from functools import partial
from Bio.Blast import NCBIXML
print = partial(print, flush=True)


# def protein_align_score(seq_a, seq_b, mode="local"):
#     aligner = Align.PairwiseAligner()
#     aligner.mode = "local"
#     score = aligner.score(a, b) / max(len(seq_a), len(seq_b))
#     return score

class PSSM(object):
    def __init__(self, fn) -> None:
        super().__init__()
        self.pssm = list()
        self.res = list()
        with open(fn) as infile:
            for nr, l in enumerate(infile):
                if nr == 1:
                    if l.startswith("Last position"):
                        self.__parse_pssm(fn)
                    else:
                        self.__parse_blossom_pssm(fn)
                    break
        self.pssm = np.concatenate(self.pssm, axis=0)
        self.res = ''.join(self.res)
    
    def __parse_pssm(self, fn):
        with open(fn) as infile:
            for nr, l in enumerate(infile):
                if nr < 3:
                    continue
                if len(l.strip()) > 1:
                    fields = re.split('\s+', l.strip())
                    self.pssm.append(np.array([int(x) for x in fields[2:22]]).reshape(1, -1))
                    self.res.append(fields[1])
                else:
                    break
    
    def __parse_blossom_pssm(self, fn):
        with open(fn) as infile:
            for nr, l in enumerate(infile):
                if nr < 1:
                    continue
                fields = re.split('\s+', l.strip())
                self.pssm.append(np.array([int(x) for x in fields[2:22]], dtype=int).reshape(1, -1))
                self.res.append(fields[1])


class SPOTDisorder(object):
    def __init__(self, fn) -> None:
        super().__init__()
        self.res = list()
        self.disorder = list()
        with open(fn) as infile:
            for l in infile:
                if l.startswith('#'):
                    continue
                _, aa, prob, _ = l.strip().split()
                self.res.append(aa)
                self.disorder.append(float(prob))
        self.res = ''.join(self.res)
        self.disorder = np.array(self.disorder, dtype=np.float16)


class DSSP(object):
    def __init__(self, fn) -> None:
        super().__init__()
        self.res= dict()
        self.ss = dict()
        self.acc = dict()
        self.__read_dssp(fn)
    
    def __read_dssp(self, fn):
        start = False
        chain = None
        with open(fn) as infile:
            for l in infile:
                if l.startswith("  #  RESIDUE AA"):
                    assert l[5:12] == "RESIDUE"
                    assert l[13:15] == "AA"
                    assert l[16:25] == "STRUCTURE"
                    assert l[35:38] == "ACC"
                    start = True
                elif start:
                    if l[13:15] == "!*":
                        chain = None
                        continue
                    elif l[13:15] == "! ":
                        continue
                    if chain is None:
                        chain = l[11]
                        self.res[chain] = list()
                        self.ss[chain] = list()
                        self.acc[chain] = list()
                    else:
                        assert chain == l[11], "{}\n{}\n{}\t{}".format(fn, l, chain, l[11])
                    aa = l[13]
                    acc = int(l[34:38].strip())
                    ss = l[16]
                    self.res[chain].append(aa)
                    self.ss[chain].append(ss)
                    self.acc[chain].append(acc)
        for c in self.acc:
            self.acc[c] = np.array(self.acc[c], dtype=np.int32)
            self.res[c] = ''.join(self.res[c])
            self.ss[c] = ''.join(self.ss[c])


class BlastMapping(object):
    def __init__(self, blast_xml: str, debug=False, version=1):
        self.query_seqs = dict()
        self.mapping = dict() # query -> sbjct -> query_range -> sbjct_range
        self.process(blast_xml, version, debug)

    # def process(self, blast_xml, query_name_fun: Callable[[str], str], sbjct_name_fun: Callable[[str], str]):
    def process(self, blast_xml, version, debug=False):
        cache = "{}.{}.pkl.gz".format(blast_xml, version)
        if not debug and os.path.exists(cache):
            print("Loading cache {} ...".format(cache), file=sys.stderr)
            self.mapping = pickle.load(gzip.open(cache, 'rb'))
        else:
            blast_records = NCBIXML.parse(gzip.open(blast_xml, 'rt'))

            for idx, record in tqdm.tqdm(enumerate(blast_records), desc="Reading {}".format(os.path.basename(blast_xml))):
                if debug and idx > 1000:
                    break
                if len(record.alignments) == 0:
                    continue
                warnings.warn("all space in query name have are removed")
                query = record.query.replace(' ', '') # e.g.: chr1_2337219_G_A|left_pad=97|right_pad=3|ref=R|alt=W

                if query in self.mapping:
                    if self.__clean_seq(record.alignments[0].hsps[0].query) != self.mapping[query]:
                        warnings.warn("Conflicted query name: {}, renamed to {}_{}".format(query, query, record.query_id))
                        query = "{}_{}".format(query, record.query_id)
                self.mapping[query] = dict()
                aligns = list()
                for a in record.alignments:
                    seq_name = a.hit_def
                    h = a.hsps[0]
                    aligns.append((-h.identities, int(a.accession), seq_name, h))
                if debug:
                    print(aligns)
                
                aligns = sorted(aligns, key=lambda x:(x[0], x[1]))
                if debug:
                    print(aligns)

                seq_name, h = aligns[0][-2:]

                self.mapping[query][seq_name] = OrderedDict()
                # h = next(iter(a.hsps))[0]
                match = h.match
                if not query and len(match.replace(' ', '').replace('+', '')) / len(match) < 0.8:
                    del self.mapping[query]
                    continue
                q_align, q_start, q_end = h.query, h.query_start, h.query_end
                s_align, s_start, s_end = h.sbjct, h.sbjct_start, h.sbjct_end
                q_start, s_start = q_start - 1, s_start - 1
                cur_query, cur_sbjct = q_start, s_start
                for i, m in enumerate(match):
                    if q_align[i] == m:
                        self.mapping[query][seq_name][cur_query] = cur_sbjct
                    if q_align[i] != '-':
                        cur_query += 1
                    if s_align[i] != '-':
                        cur_sbjct += 1
            if not debug:
                pickle.dump(self.mapping, gzip.open(cache, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    
    def __clean_seq(self, s):
        s = s.strip().replace('-', '').replace('-', '')
        return s
    

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('vcf')
    p.add_argument('--mut-protein', default="./final_train_070821.fixed.protein.mark_mut.fa.gz", required=True)
    p.add_argument('--blast', default="./refseq2uniprot.blastp.xml.gz", help="blast results", required=True)
    p.add_argument('-o', required=True)

    # p.add_argument('-vcf', default="./final_train_070821.fixed.marked.vcf")
    # p.add_argument('--mut-protein', default="./final_train_070821.fixed.protein.mark_mut.fa.gz")
    # p.add_argument('--blast', default="./refseq2uniprot.blastp.xml.gz", help="blast results")
    #p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)

    # line2name = dict()
    # with open(args.vcf) as infile:
    #     for nr, l in enumerate(infile):
    #         name = l.split('\t')[2].split('|')[0]
    #         line2name["line{}".format(nr + 1)] = name
    
    # index_mapping, pairwise_identity = load_seq_mapping(args.blast)

    # aa2idx = dict()
    # for i, aa in enumerate("ARNDCQEGHILKMFPSTWYV"):
    #     aa2idx[aa] = i

    # pro_pssms = dict()
    # pro_disorder = dict()

    # total_size = 50
    # sample_pssm = dict()
    # with gzip.open(args.mut_protein, 'rt') as infile:
    #     for l in tqdm.tqdm(infile):
    #         if l.startswith('>'):
    #             name = l.strip().lstrip('>')
    #         else:
    #             tx_id = name.split()[1]
    #             sample_id = line2name[name.split(' ')[0]]

    #             left_ref, alt_right = l.strip().split('/')
    #             left, ref = left_ref.split('[')
    #             alt, right = alt_right.split(']')
    #             ref, alt = ref.strip('*'), alt.strip('*')
    #             ref_seq = ''.join((left, ref, right)).strip("*")
    #             alt_seq = ''.join((left, alt, right)).strip("*")

    #             if tx_id not in index_mapping:
    #                 continue

    #             try:
    #                 pro_id = next(iter(index_mapping[tx_id].keys()))
    #             except:
    #                 warnings.warn("missing uniprot sequence: {}".format(tx_id))
    #                 continue
    #             if pro_id not in pro_pssms:
    #                 pssm_fn = "/home/chenken/db/PSSM-HHM-SPOT_DISORDER/{}/{}.pssm".format(pro_id[-2:], pro_id)
    #                 disorder_fn = "/home/chenken/db/PSSM-HHM-SPOT_DISORDER/{}/{}.spotds".format(pro_id[-2:], pro_id)
    #                 if not os.path.exists(pssm_fn) or not os.path.exists(disorder_fn):
    #                     warnings.warn("{} missing profile: {}".format(tx_id, pro_id))
    #                     continue
    #                 pro_pssms[pro_id] = load_pssm(pssm_fn)
    #                 pro_disorder[pro_id] = load_disorder(disorder_fn)
    #             pro_pssm = pro_pssms[pro_id]
    #             pro_diso = pro_disorder[pro_id]
    #             if len(left) > len(pro_pssm):
    #                 warnings.warn("{} exceeds the length of {}({})".format(len(left), pro_id, len(pro_pssm)))

    #             mut_size = min(max(len(ref), len(alt)), total_size)
    #             left_pad = min(max(0, (total_size - mut_size) // 2), len(left))
    #             assert left_pad >= 0

    #             ref_pssm = list()
    #             alt_mark = list()
    #             ref_diso = list()

    #             idx_mapping = index_mapping[tx_id][pro_id]

    #             skip = False
    #             for i in range(left_pad):
    #                 idx = len(left) - left_pad + i
    #                 if idx not in idx_mapping:
    #                     skip = True
    #                     break
    #                 idx = idx_mapping[idx]
    #                 ref_pssm.append(pro_pssm[idx])
    #                 ref_diso.append(pro_diso[idx])
    #                 alt_mark.append(np.array([0, 0, 0], dtype=int).reshape(-1, 1))

    #             if not skip:
    #                 for i in range(min(mut_size, total_size - left_pad)):
    #                     idx = len(left) + i
    #                     if i < len(ref) and i < len(alt):
    #                         if idx not in idx_mapping:
    #                             skip = True
    #                             break
    #                         idx = idx_mapping[idx]
    #                         try:
    #                             pssm_ar = pro_pssm[idx]
    #                         except:
    #                             print(name, pro_id, idx, file=sys.stderr)
    #                             raise IndexError("{}".format((i, len(left), left_pad, idx, len(pro_pssm))))
    #                         ref_pssm.append(pssm_ar)
    #                         ref_diso.append(pro_diso[idx])
    #                         nr, na = ref[i], alt[i]
    #                         # try:
    #                         change = pssm_ar[aa2idx[na]] - pssm_ar[aa2idx[nr]]
    #                         # except:
    #                         #     print("{}\n{}\n{}\n{}/{}".format(name, ref_seq, alt_seq, nr, na), file=sys.stderr)
    #                         #     raise KeyError()
    #                         alt_mark.append(np.array([0, 0, change], dtype=int).reshape(-1, 1))
    #                     elif i < len(ref):
    #                         if idx not in idx_mapping:
    #                             skip = True
    #                             break
    #                         idx = idx_mapping[idx]
    #                         pssm_ar = pro_pssm[idx]
    #                         ref_pssm.append(pssm_ar)
    #                         ref_diso.append(pro_diso[idx])
    #                         alt_mark.append(np.array([1, 0, 0], dtype=int).reshape(-1, 1))
    #                     elif i < len(alt):
    #                         ref_pssm.append(np.zeros(20).reshape(-1, 1))
    #                         ref_diso.append(0)
    #                         alt_mark.append(np.array([0, 1, 0], dtype=int).reshape(-1, 1))

    #             if not skip:
    #                 for i in range(max(total_size - mut_size - left_pad, 0)):
    #                     idx = i + len(left) + mut_size
    #                     if idx not in idx_mapping:
    #                         skip = True
    #                         break
    #                     idx = idx_mapping[idx]
    #                     ref_pssm.append(pro_pssm[idx])
    #                     ref_diso.append(pro_diso[idx])
    #                     alt_mark.append(np.array([0, 0, 0]).reshape(-1, 1))

    #             if not skip:
    #                 try:
    #                     ref_pssm = np.concatenate(ref_pssm, axis=1)
    #                     ref_diso = np.array(ref_diso, dtype=np.float16).reshape(1, -1)
    #                     alt_mark = np.concatenate(alt_mark, axis=1)
    #                 except:
    #                     print(ref_pssm, file=sys.stderr)
    #                     print(alt_mark, file=sys.stderr)
    #                     raise ValueError()

    #                 sample_pssm[sample_id] = np.concatenate((ref_diso, ref_pssm, alt_mark), axis=0)
    #                 # print("{}\t{}".format(sample_id, sample_pssm[sample_id].shape))

    # torch.save(sample_pssm, args.o)


