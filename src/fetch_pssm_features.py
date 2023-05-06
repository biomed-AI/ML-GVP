#!/usr/bin/env python3

import argparse, os, sys, time, tqdm, re, pickle

from numpy import dtype, lib
from utils import copen
import warnings, json, gzip, pickle
import logging
logger = logging.getLogger(__name__)
from multiprocessing import Pool, process
# warning = logging.warning

import numpy as np
from typing import Any, Dict, List, Union
from functools import partial
print = partial(print, flush=True)
from collections import defaultdict, OrderedDict
import libprotmapping

AA2IDX = {aa: i for i, aa in enumerate("ARNDCQEGHILKMFPSTWYV")}

AA_ONEHOT = dict()
for i, aa in enumerate("ARNDCQEGHILKMFPSTWYV"):
    ar = [0] * 20
    ar[i] += 10
    AA_ONEHOT[aa] = np.array(ar).reshape(-1, 1)

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--input', 
            required=True,
        # default=os.path.join(DATA_DIR, "final_test_092921.raw.avinput.coding_change.marked.substr.fa"),
        help="substr protein fasta, e.g.: final_test_092921.protein.substr.fa")
    p.add_argument('-b', "--blast-map", 
        # default=os.path.join(DATA_DIR, "final_test_092921.coding_change.marked.substr.uniprot.xml.gz"), 
        required=True,
        help="Blastp output in xml format")
    p.add_argument("-db", "--path-to-pssm-spot", required=True)
    p.add_argument('--seed', type=int, default=2020)
    return p


def prepare_mark(samples):
    marks = list()
    # for i in range(len(left_pads)):
    for name, sample_id, left_pad, right_pad, ref, alt, ref_prot_seq in samples:
        mark = np.concatenate((
            np.arange(1, left_pad + 1)[::-1],
            np.zeros(max(len(ref), len(alt))),
            np.arange(1, right_pad + 1)
        ), axis=0)
        marks.append(np.arcsinh(mark / 100).astype(np.float16).reshape(1, -1))
    return marks


def prepare_dssp_disorder(samples, blast_map: libprotmapping.BlastMapping, pssm_spot_db):
    pssm_cache, disorder_cache = dict(), dict()
    for name, sample_id, left_pad, right_pad, ref, alt, ref_prot_seq in tqdm.tqdm(samples, desc="samples"):
        name = name.replace(' ', '')
        if name in blast_map.mapping:
            ref_id = iter(blast_map.mapping[name]).__next__()
        else:
            print("not found: {}".format(name))
            continue
        # if ref_id.startswith("NM_"):
        #     ref_name = ref_id
        #     pssm_fn = "/home/chenken/Documents/CAGI6-Sherloc/data/cagi6/refseqs/{}.pssm".format(ref_name)
        #     disorder_fn = "/home/chenken/Documents/CAGI6-Sherloc/data/cagi6/refseqs/{}.spotds".format(ref_name)
        # else:
        ref_name = ref_id.split('|')[1]
        pssm_fn ="{}/{}/{}.pssm".format(pssm_spot_db, ref_name[-2:], ref_name) 
        disorder_fn = "{}/{}/{}.spotds".format(pssm_spot_db, ref_name[-2:], ref_name)

        if ref_name in pssm_cache:
            pssm_profile = pssm_cache[ref_name]
            disorder_profile = disorder_cache[ref_name]
        elif os.path.exists(pssm_fn):
            pssm_profile = libprotmapping.PSSM(pssm_fn)
            disorder_profile = libprotmapping.SPOTDisorder(disorder_fn)
            pssm_cache[ref_name] = pssm_profile
            disorder_cache[ref_name] = disorder_profile
        else:
            pssm_profile = None
            disorder_fn = None

        missing, keep = 0, True 
        cont_disorder, cont_pssm, cont_pssm_change = list(), list(), list()
        point_disorder, point_pssm, point_pssm_change = list(), list(), list()
        for i in range(left_pad):
            idx = i
            if idx in blast_map.mapping[name][ref_id]:
                j = blast_map.mapping[name][ref_id][idx]
                if j >= len(disorder_profile.res) or idx >= len(ref_prot_seq) or disorder_profile.res[j] != ref_prot_seq[idx]:
                    print("({},{})\n{}\n{}\n{}".format(
                        disorder_fn, pssm_fn, 
                        ref_prot_seq, 
                        disorder_profile.res, 
                        blast_map.mapping[name][ref_id]
                    ))
                    raise RuntimeError()
                cont_disorder.append(disorder_profile.disorder[j])
                point_disorder.append(disorder_profile.disorder[j])
                cont_pssm.append(pssm_profile.pssm[j].reshape(-1, 1))
            else:
                missing += 1
                # warnings.warn("Missing residue {}({}) in {} -> {}".format(ref_prot_seq[i], i, sample_id, ref_name))
                cont_disorder.append(cont_disorder[-1] if i > 0 else 0.5)
                point_disorder.append(cont_disorder[-1] if i > 0 else 0.5)
                cont_pssm.append(AA_ONEHOT[ref_prot_seq[i]])
            cont_pssm_change.append(np.array([0, 0, 0]).reshape(-1, 1))
        for i in range(min(len(ref), len(alt))):
            idx = i + left_pad
            if idx in blast_map.mapping[name][ref_id]:
                j = blast_map.mapping[name][ref_id][idx]
                cont_disorder.append(disorder_profile.disorder[j])
                point_disorder.append(disorder_profile.disorder[j])
                cont_pssm.append(pssm_profile.pssm[j].reshape(-1, 1))
                cont_pssm_change.append(np.array([
                    0, 
                    0, 
                    pssm_profile.pssm[j][AA2IDX[alt[i]]] - pssm_profile.pssm[j][AA2IDX[ref[i]]]
                ]).reshape(-1, 1))
                point_pssm.append(pssm_profile.pssm[j].reshape(-1, 1))
                point_pssm_change.append(abs(pssm_profile.pssm[j][AA2IDX[alt[i]]] - pssm_profile.pssm[j][AA2IDX[ref[i]]]))
            else:
                point_disorder.append(0.5)
                missing += 1
                keep = False
                # break
        if len(ref) > len(alt):
            for i in range(len(ref) - len(alt)):
                idx = i + left_pad + len(alt)
                point_pssm_change.append(30)
                if idx in blast_map.mapping[name][ref_id]: 
                    j = blast_map.mapping[name][ref_id][idx]
                    cont_disorder.append(disorder_profile.disorder[j])
                    point_disorder.append(disorder_profile.disorder[j])
                    cont_pssm.append(pssm_profile.pssm[j].reshape(-1, 1))
                    point_pssm.append(pssm_profile.pssm[j].reshape(-1, 1))
                    cont_pssm_change.append(np.array([1, 0, 0]).reshape(-1, 1)) # deletion
                else:
                    missing += 1
                    cont_disorder.append(cont_disorder[-1] if idx > 0 else 0.5)
                    point_disorder.append(cont_disorder[-1] if idx > 0 else 0.5)
                    cont_pssm.append(AA_ONEHOT[ref_prot_seq[idx]])
                    cont_pssm_change.append(np.array([0, 0, 0]).reshape(-1, 1))
        elif len(ref) < len(alt):
            for i in range(len(alt) - len(ref)):
                point_pssm_change.append(30)
                cont_disorder.append(0)
                cont_pssm.append(np.zeros((20, 1)))
                cont_pssm_change.append(np.array([0, 1, 0]).reshape(-1, 1)) # insertion 

        for i in range(right_pad):
            idx = i + left_pad + len(ref)
            if idx in blast_map.mapping[name][ref_id]:
                j = blast_map.mapping[name][ref_id][idx]
                cont_disorder.append(disorder_profile.disorder[j])
                point_disorder.append(disorder_profile.disorder[j])
                cont_pssm.append(pssm_profile.pssm[j].reshape(-1, 1))
            else:
                missing += 1
                cont_disorder.append(cont_disorder[-1] if len(cont_disorder) > 0 else 0.5)
                point_disorder.append(cont_disorder[-1] if len(cont_disorder) > 0 else 0.5)
                cont_pssm.append(AA_ONEHOT[ref_prot_seq[idx]])
            cont_pssm_change.append(np.array([0, 0, 0]).reshape(-1, 1))

        left, right = left_pad, left_pad + len(ref)
        assert len(point_disorder) == left_pad + len(ref) + max(0, right_pad), "{}".format((name, len(point_disorder), left_pad, right_pad, ref, alt))
        point_disorder = (
                point_disorder[left:right] if right > left else point_disorder[max(left - 1, 0):right+1], 
                point_disorder[max(0, left-3):left + len(ref) + 3], 
                point_disorder[max(0, left-10): left + len(ref) + 10]
            )
        if len(point_pssm) == 0:
            warnings.warn(name)
            point_pssm = np.array([np.zeros(20)])
        else:
            point_pssm = np.concatenate(point_pssm, axis=1).T
        prot_features[sample_id] = {
                "disorder": point_disorder,
                "pssm": point_pssm.astype(np.float16),
                "pssm_change": point_pssm_change
                }

        # if not keep:
        #     logging.warning("Remove: {}".format(sample_id))
        # elif missing / len(ref_prot_seq) > 0.5:
        #     logging.warning("Remove: {} (missing rate {:.3f}, {}, {})".format(sample_id, missing / len(ref_prot_seq), sample_id, ref_name))
        # else:
        #     cont_disorder = np.array(disorder, dtype=np.float16).reshape(1, -1)
        #     cont_pssm = (np.concatenate(pssm, axis=1) / 10).astype(np.float16)
        #     cont_pssm_change = (np.concatenate(pssm_change, axis=1) / 10).astype(np.float16)
        #     prot_features[sample_id] = {
        #         "cont_disorder": cont_disorder,
        #         "cont_pssm": cont_pssm,
        #         "cont_pssm_change": cont_pssm_change
        #     }
    return prot_features


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    np.random.seed(args.seed)


    samples = list()
    prot_features = dict()
    with copen(args.input) as infile:
        for l in tqdm.tqdm(infile, desc="Reading {} ...".format(os.path.basename(args.input))):
            if l.startswith('>'):
                # name, left_pad, right_pad, ref, alt = l.strip().split('>').split('|')
                name = l.strip().strip('>')
            else:
                ref_prot_seq = l.strip().strip('*')
                sample_id, _, _, left_right_pad, ref_alt = name.split('|')
                left_pad, right_pad = left_right_pad.split('=')[-1].split('/')
                left_pad, right_pad = int(left_pad), int(right_pad)
                ref, alt = ref_alt.split('=')[-1].split('/')
                samples.append([name, sample_id, left_pad, right_pad, ref, alt, ref_prot_seq])

                # prot_features[name] = {
                #     "metainfo": (
                #         sample_id, 
                #         left_pad, right_pad, 
                #         ref, alt, 
                #         ref_prot_seq
                #     ),

                # }


                # pdb_ref = iter(blast_pdb[sample_id]).__next__()
                # alphafold_ref = iter(blast_alphafold[sample_id]).__next__()

                # ## secondary structure
                # ss = list()
                # asa = list()
    # n_batch = 16
    # batch_size = np.ceil(len(samples) / n_batch).astype(int)
    # job_list = list()
    # for i in range(n_batch):
    #     job_list.append([
    #         samples[i * batch_size:(i + 1) * batch_size], 
    #         blast_uniprot, 
    #         i
    #     ])
    # with Pool(processes=n_batch) as pool:
    #     res = pool.starmap(prepare_dssp_disorder, job_list)
    # logging.info("Merging ... {}".format(time.asctime()))
    # prot_features = list()
    # for i in range(n_batch - 1, -1, -1):
    #     prot_features.extend(res[i])
    #     del res[i]
    # logging.info("Saving ... {}".format(time.asctime()))

    # blast_refseq = libprotmapping.BlastMapping(args.blast_to_refseq)
    # blast_uniprot = libprotmapping.BlastMapping(args.blast_to_uniprot)
    blast_map = libprotmapping.BlastMapping(args.blast_map)

    prot_features = prepare_dssp_disorder(samples, blast_map, pssm_spot_db=args.path_to_pssm_spot)
    pickle.dump(
            prot_features, 
            gzip.open("{}.pssm-disorder.pkl.gz".format(args.input), 'wb'), 
            protocol=pickle.HIGHEST_PROTOCOL
        )

# def prepare_ss(samples, blast_pdb: libprotmapping.BlastMapping, blast_alphafold: libprotmapping.BlastMapping):
