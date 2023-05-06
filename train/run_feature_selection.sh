#!/bin/bash

CUDA_VISIBLE_DEVICES=3 ./feature_selection.py \
    --data ./features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
    -o ./feature_configs/feature_selection.nonpro-other.json \
    --device 0 \
    -f ./feature_configs/no_protein_features.json > ../log/feature_selection.nonpro-other.log