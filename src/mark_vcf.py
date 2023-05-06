#!/usr/bin/env python3

import argparse
from collections import defaultdict
from utils import copen
import json
from functools import partial
print = partial(print, flush=True)


def strip_ENS_version(ensembl_id: str) -> str:
    suffix = ""
    if ensembl_id.endswith("_PAR_Y"):
        suffix = "_PAR_Y"
    return "{}{}".format(ensembl_id.split('.')[0], suffix)

pro_change_re = {
    "missense": "^p.([A-Z])(\d+)([A-Z])"
}

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('vcf')
    p.add_argument('-f', '--variant_function', required=True)
    p.add_argument('-e', "--exonic_variant_function", required=True)
    p.add_argument('-m', "--tx2uniprot", required=True, help="transcript -> uniprot mapping file")
    p.add_argument('-tx', required=True, help="tss bed")
    return p

def load_tx_info(fn):
    tx_info = dict()
    with copen(fn) as infile:
        for l in infile:
            gene_id, gene_name, gene_type, tx_id, _, strand = l.split('\t')[3].split('|')
            tx_id = strip_ENS_version(tx_id)
            tx_info[tx_id] = (gene_id, gene_name, gene_type, strand)
    return tx_info


if __name__ == "__main__":
    p = get_args()
    args = p.parse_args()
    #np.random.seed(args.seed)

    tx2uniprot = json.load(open(args.tx2uniprot))

    exonic_function_info = dict()
    tx_cnt = defaultdict(int)
    with open(args.exonic_variant_function) as infile:
        for l in infile:
            nline, c_change, mut_info = l.strip().split('\t')[0:3]
            nline = int(nline.replace("line", '')) - 1
            mut_info = mut_info.split(',')
            assert len(mut_info) == 2 and mut_info[1] == "", "{}".format(mut_info)
            mut_info = mut_info[0].split(':')
            tx_id = mut_info[1]
            consequence = ':'.join(mut_info[2:])
            tx_cnt[tx_id] += 1
            exonic_function_info[nline] = (tx_id, consequence)

    tx_info = load_tx_info(args.tx)
    variant_annotations = list()
    with open(args.variant_function) as infile:
        for nr, l in enumerate(infile):
            region, txs = l.strip().split('\t')[0:2]
            if nr in exonic_function_info:
                tx_id, consequence = exonic_function_info[nr]
            else:
                consequence = "others"
                tx_list = list()
                for tx in [t .split('(')[0] for t in txs.replace(';', ',').split(',')]:
                    if tx not in tx2uniprot:
                        tx_list.append((tx, 0, tx_cnt[tx]))
                    else:
                        tx_list.append((tx, tx2uniprot[tx][1], tx_cnt[tx]))
                tx_id = sorted(tx_list, key=lambda ar:(ar[1], ar[2]), reverse=True)[0][0]
            try:
                gene_id, gene_name, gene_type, strand = tx_info[tx_id]
                variant_annotations.append((
                    region, tx_id, consequence, gene_id, gene_name, gene_type, strand
                ))
            except:
                gene_id, gene_name, gene_type, strand = '.', '.', '.', '+'
                variant_annotations.append((
                    region, tx_id, consequence, gene_id, gene_name, gene_type, strand
                ))


    with open(args.vcf) as infile:
        nr = 0
        for _, l in enumerate(infile):
            if l.startswith("#"):
                continue
            chrom, pos, mut_id, ref, alt = l.strip().split('\t')[0:5]
            try:
                mut_id = "{}@@{}".format(mut_id, '@@'.join(variant_annotations[nr]))
                print('\t'.join([chrom, pos, mut_id, ref, alt]))
            except:
                pass
            nr += 1
        pass


    # updated_tx = {
    #     "NM_001257360": "NM_001368809",
    #     "NM_032374": "NM_001370595",
    #     "NM_000060": "NM_001370658"
    # }

    # refseq2uniprot = json.load(open(args.refseq2uniprot))
    # uniprot2pdb = json.load(open(args.uniprot2pdb))
    # transcripts = load_transcripts(fn=args.db)

    # variant_function = list()
    # with open(args.variant_function) as infile:
    #     for nr, l in enumerate(infile):
    #         variant_function.append(l.split('\t')[0:2])
    # exonic_variant_function = [("non-coding", "non-coding") for _ in variant_function]
    # with open(args.exonic_variant_function) as infile:
    #     for l in infile:
    #         line, mut_type, mut_info = l.strip().split('\t')[0:3]
    #         line = int(line.replace("line", "")) - 1
    #         exonic_variant_function[line] = (mut_type.replace(' ', '_'), mut_info)

    # nlines = list()
    # with open(args.vcf) as infile:
    #     for l in infile:
    #         nlines.append(l)
    # assert len(variant_function) == len(nlines)
    
    # re_tx_id = re.compile(r'N[MR]_\d+')

    # with open(args.vcf) as infile:
    #     for nr, l in enumerate(infile):
    #         chrom, position, name, ref, alt = l.strip().split('\t')[0:5]
    #         
    #         raw_tx = name.split('|')[1].split('.')[0]
    #         if raw_tx in updated_tx:
    #             raw_tx = updated_tx[raw_tx]
    #         region, anno_txs = variant_function[nr]

    #         txs = set(re_tx_id.findall(anno_txs))
    #         # if raw_tx not in txs:
    #         #     logging.warning("{}\t{}".format(raw_tx, anno_txs))
    #         if "{}@{}".format(raw_tx, chrom) in transcripts:
    #             tx_id = raw_tx
    #         else:
    #             if len(txs) > 0:
    #                 max_len = 0
    #                 for tx in sorted(txs):
    #                     if transcripts["{}@{}".format(tx, chrom)].tx_length > max_len:
    #                         tx_id = tx
    #                         max_len = transcripts["{}@{}".format(tx, chrom)].tx_length
    #                 logging.warning("multiple\t{}\t{}\t{} -> {}".format(raw_tx, region, anno_txs, tx_id))
    #         # try:
    #         #     _ = tx_id
    #         # except:
    #         #     raise RuntimeError("{}".format(raw_tx))
    #         strand = transcripts["{}@{}".format(tx_id, chrom)].strand

    #         pdb = set()
    #         if tx_id in refseq2uniprot:
    #             uniprot_ids = ','.join(sorted(refseq2uniprot[tx_id]))
    #             for u in uniprot_ids.split(','):
    #                 if u in uniprot2pdb:
    #                     pdb = pdb.union(set(uniprot2pdb[u]))
    #         else:
    #             uniprot_ids = "no_uniprot_found"
    #         if len(pdb) == 0:
    #             pdb = ["no_pdb_found"]
    #         pdb = ','.join(sorted(list(pdb)))
    #         # print("{}\t{}\t{}|{}|{}\t{}\t{}".format(chrom, position, name, tx_id, strand, ref, alt))
    #         print("{}\t{}\t{}|{}|{}|{}|{}|{}|{}|{}\t{}\t{}".format(chrom, position, name, tx_id, strand, variant_function[nr][0], exonic_variant_function[nr][0], exonic_variant_function[nr][1], uniprot_ids, pdb, ref, alt))
    #         del tx_id, raw_tx

