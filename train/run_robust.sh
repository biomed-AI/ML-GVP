#!/bin/bash


if [ $# -lt 1 ]; then
    echo "usage: $0 mode(null, rfecv)"
    exit 1
fi

if [ ! $1 ]
then
    option=$1
    features_protein_snp="all_features.yml"
    features_protein_other="all_features.yml"
    features_nonpro_snp="no_protein_features.yml"
    features_nonpro_other="no_protein_features.yml"
elif [ "$1" == "rfecv" ]
then
    option=_$1
    features_protein_snp="feature_selection.protein-snp.yml"
    features_protein_other="feature_selection.protein-other.yml"
    features_nonpro_snp="feature_selection.nonpro-snp.yml"
    features_nonpro_other="feature_selection.nonpro-other.yml"
else
    echo "mode: null or rfecv"
    exit 1
fi
WORKING_DIR=`dirname $0`

# coding
export CUDA_VISIBLE_DEVICES=0
./robust.py \
    -c $WORKING_DIR/model_configs/bayes_opt_protein-snp${option}.json \
    --training_data $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
    --test_data $WORKING_DIR/../final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical.protein_snp.tsv \
    -p robust_protein_snp${option} \
    -f $WORKING_DIR/feature_configs/$features_protein_snp


export CUDA_VISIBLE_DEVICES=0
./robust.py \
    -c $WORKING_DIR/model_configs/bayes_opt_protein-other${option}.json \
    --training_data $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
    --test_data $WORKING_DIR/../final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical.protein_other.tsv \
    -p robust_protein_other${option} \
    -f $WORKING_DIR/feature_configs/$features_protein_other


export CUDA_VISIBLE_DEVICES=2
./robust.py \
    -c $WORKING_DIR/model_configs/bayes_opt_nonpro-snp${option}.json \
    --training_data $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
    --test_data $WORKING_DIR/../final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical.nonpro_snp.tsv \
    -p robust_nonpro_snp${option} \
    -f $WORKING_DIR/feature_configs/$features_nonpro_snp


export CUDA_VISIBLE_DEVICES=2
./robust.py \
    -c $WORKING_DIR/model_configs/bayes_opt_nonpro-other${option}.json \
    --training_data $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
    --test_data $WORKING_DIR/../final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical.nonpro_other.tsv \
    -p robust_nonpro_other${option} \
    -f $WORKING_DIR/feature_configs/$features_nonpro_other


# export CUDA_VISIBLE_DEVICES=0
# ./robust.py \
#     -c $WORKING_DIR/model_configs/bayes_opt_protein-snp${option}.json \
#     --training_data $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
#     --test_data $WORKING_DIR/../final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical.protein_snp.tsv \
#     -p robust_protein_snp${option} \
#     -f $WORKING_DIR/feature_configs/$features_protein_snp &> ../log/robust_protein_snp${option}.log &


# export CUDA_VISIBLE_DEVICES=0
# ./robust.py \
#     -c $WORKING_DIR/model_configs/bayes_opt_protein-other${option}.json \
#     --training_data $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
#     --test_data $WORKING_DIR/../final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical.protein_other.tsv \
#     -p robust_protein_other${option} \
#     -f $WORKING_DIR/feature_configs/$features_protein_other & > ../log/robust_protein_other${option}.log &


# export CUDA_VISIBLE_DEVICES=2
# ./robust.py \
#     -c $WORKING_DIR/model_configs/bayes_opt_nonpro-snp${option}.json \
#     --training_data $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
#     --test_data $WORKING_DIR/../final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical.nonpro_snp.tsv \
#     -p robust_nonpro_snp${option} \
#     -f $WORKING_DIR/feature_configs/$features_nonpro_snp & > ../log/robust_nonpro_snp${option}.log &


# export CUDA_VISIBLE_DEVICES=2
# ./robust.py \
#     -c $WORKING_DIR/model_configs/bayes_opt_nonpro-other${option}.json \
#     --training_data $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
#     --test_data $WORKING_DIR/../final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical.nonpro_other.tsv \
#     -p robust_nonpro_other${option} \
#     -f $WORKING_DIR/feature_configs/$features_nonpro_other & > ../log/robust_nonpro_other${option}.log

# wait

# ../src/merge_4_model.py \
#     --ps $WORKING_DIR/predictions/robust_protein_snp${option}.prediction.txt \
#     --po $WORKING_DIR/predictions/robust_protein_other${option}.prediction.txt \
#     --ns $WORKING_DIR/predictions/robust_nonpro_snp${option}.prediction.txt \
#     --no $WORKING_DIR/predictions/robust_nonpro_other${option}.prediction.txt  > $WORKING_DIR/predictions/train${option}.merged.tsv

# wait

# export CUDA_VISIBLE_DEVICES=3
# ./linearregression.py \
#     -d $WORKING_DIR/predictions/train${option}.merged.tsv \
#     -p final${option} \
#     -f $WORKING_DIR/feature_configs/four_group_features.yml & > ../log/final${option}.log
