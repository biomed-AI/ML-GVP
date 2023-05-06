#!/bin/bash



## coding snp
export CUDA_VISIBLE_DEVICES=3
./grid_cv_improved.py \
    --device 3 \
    -d train.protein-snp.tsv \
    -f ./all_features.json \
    -o bayes_opt_protein_snp.json &> bayes_opt_protein_snp.json.log &

## coding snp
export CUDA_VISIBLE_DEVICES=3
./grid_cv_improved.py \
    --device 3 \
    -d train.protein-other.tsv \
    -f ./all_features.json \
    -o bayes_opt_protein_other.json &> bayes_opt_protein_other.json.log &


## coding snp
export CUDA_VISIBLE_DEVICES=3
./grid_cv_improved.py \
    --device 3 \
    -d train.nonpro-snp.tsv \
    -f ./no_protein_features.json \
    -o bayes_opt_nonpro_snp.json &> bayes_opt_nonpro_snp.json.log &

## coding snp
export CUDA_VISIBLE_DEVICES=3
./grid_cv_improved.py \
    --device 3 \
    -d train.nonpro-other.tsv \
    -f ./no_protein_features.json \
    -o bayes_opt_nonpro_other.json &> bayes_opt_nonpro_other.json.log &




wait
 
#key	label	group	chrom	position	ref	alt	region	splicing	coding
