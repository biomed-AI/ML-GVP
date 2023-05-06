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


#if [ "$option" == "_rfecv" ];
#then
#    features_protein_snp="feature_selection.protein-snp.yml"
#    features_protein_other="feature_selection.protein-other.yml"
#    features_nonpro_snp="feature_selection.nonpro-snp.yml"
#    features_nonpro_other="feature_selection.nonpro-other.yml"
#else
#    features_protein_snp="all_features.yml"
#    features_protein_other="all_features.yml"
#    features_nonpro_snp="no_protein_features.yml"
#    features_nonpro_other="no_protein_features.yml"
#fi

# coding


export CUDA_VISIBLE_DEVICES=3
./cv.py \
    -c $WORKING_DIR/model_configs/bayes_opt_protein-snp${option}.json \
    -d $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
    -p cv_protein_snp${option} \
    -f $WORKING_DIR/feature_configs/$features_protein_snp &> ../log/cv_protein_snp${option}.log &


export CUDA_VISIBLE_DEVICES=3
./cv.py \
    -c $WORKING_DIR/model_configs/bayes_opt_protein-other${option}.json \
    -d $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
    -p cv_protein_other${option} \
    -f $WORKING_DIR/feature_configs/$features_protein_other & > ../log/cv_protein_other${option}.log &


export CUDA_VISIBLE_DEVICES=3
./cv.py \
    -c $WORKING_DIR/model_configs/bayes_opt_nonpro-snp${option}.json \
    -d $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
    -p cv_nonpro_snp${option} \
    -f $WORKING_DIR/feature_configs/$features_nonpro_snp & > ../log/cv_nonpro_snp${option}.log &


export CUDA_VISIBLE_DEVICES=3
./cv.py \
    -c $WORKING_DIR/model_configs/bayes_opt_nonpro-other${option}.json \
    -d $WORKING_DIR/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
    -p cv_nonpro_other${option} \
    -f $WORKING_DIR/feature_configs/$features_nonpro_other & > ../log/cv_nonpro_other${option}.log

wait

../src/merge_4_model.py \
    --ps $WORKING_DIR/predictions/cv_protein_snp${option}.prediction.txt \
    --po $WORKING_DIR/predictions/cv_protein_other${option}.prediction.txt \
    --ns $WORKING_DIR/predictions/cv_nonpro_snp${option}.prediction.txt \
    --no $WORKING_DIR/predictions/cv_nonpro_other${option}.prediction.txt  > $WORKING_DIR/predictions/train${option}.merged.tsv

wait

export CUDA_VISIBLE_DEVICES=3
./linearregression.py \
    -d $WORKING_DIR/predictions/train${option}.merged.tsv \
    -p final${option} \
    -f $WORKING_DIR/feature_configs/four_group_features.yml & > ../log/final${option}.log
