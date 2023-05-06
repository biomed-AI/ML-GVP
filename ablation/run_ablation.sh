#!/bin/bash

if [ $# -lt 1 ]; then
    echo "usage: $0 prefix"
    exit 1
fi

files=`ls ./features_config -tr|grep feature_selection.nonpro-other`
WORKDING_DIR=`dirname $0`
WORKDING_DIR=$WORKDING_DIR/../

pre=$1
# pre=`basename $input`
# pre=${pre%.gz}
# pre=${pre%_noncoding.tsv} #prefix=final_test_092921.raw.std.avinput.features-ALL-clinvar_20220730_critical


export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
    -c $WORKDING_DIR/train/model_configs/bayes_opt_protein-snp_rfecv.json \
    -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
    -p cv_protein_snp \
    -f features_config/feature_selection.protein-snp.yml &> log/cv_protein_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
    -c $WORKDING_DIR/train/model_configs/bayes_opt_protein-other_rfecv.json \
    -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
    -p cv_protein_other \
    -f features_config/feature_selection.protein-other.yml &> log/cv_protein_other.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
    -c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-snp_rfecv.json \
    -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
    -p cv_nonpro_snp \
    -f features_config/feature_selection.nonpro-snp.yml &> log/cv_nonpro_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
    -c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-other_rfecv.json \
    -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
    -p cv_nonpro_other \
    -f features_config/feature_selection.nonpro-other.yml & > log/cv_nonpro_other.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/cv_protein_snp.prediction.txt \
    --po predictions/cv_protein_other.prediction.txt \
    --ns predictions/cv_nonpro_snp.prediction.txt \
    --no predictions/cv_nonpro_other.prediction.txt  > predictions/train.merged.tsv
wait
export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/linearregression.py \
    -d predictions/train.merged.tsv \
    -p final_train \
    -f $WORKDING_DIR/train/feature_configs/four_group_features.yml > log/final.log
    
for type in "noncoding" "coding" "missense" "synonymous";
do
    echo " - "${type}
    prefix=${pre}_${type}
    $WORKDING_DIR/src/split_data.py $prefix.tsv -p $prefix
    $WORKDING_DIR/src/predict.py \
        -i ${prefix}.protein_snp.tsv \
        -m $WORKDING_DIR/models/cv_protein_snp.model.pkl \
        -p predictions/${prefix}.protein_snp &> log/${prefix}.protein_snp.log &

    $WORKDING_DIR/src/predict.py \
        -i ${prefix}.protein_other.tsv \
        -m $WORKDING_DIR/models/cv_protein_other.model.pkl \
        -p predictions/${prefix}.protein_other &> log/${prefix}.protein_other.log &

    $WORKDING_DIR/src/predict.py \
        -i ${prefix}.nonpro_snp.tsv \
        -m $WORKDING_DIR/models/cv_nonpro_snp.model.pkl \
        -p predictions/${prefix}.nonpro_snp &> log/${prefix}.nonpro_snp.log &

    $WORKDING_DIR/src/predict.py \
        -i ${prefix}.nonpro_other.tsv \
        -m $WORKDING_DIR/models/cv_nonpro_other.model.pkl \
        -p predictions/${prefix}.nonpro_other &> log/${prefix}.nonpro_other.log &
    wait

    $WORKDING_DIR/src/merge_4_model.py \
        --ps predictions/${prefix}.protein_snp.prediction.txt \
        --po predictions/${prefix}.protein_other.prediction.txt \
        --ns predictions/${prefix}.nonpro_snp.prediction.txt \
        --no predictions/${prefix}.nonpro_other.prediction.txt > predictions/test_${type}.merged.tsv

    $WORKDING_DIR/src/predict_lr.py \
        -i predictions/test_${type}.merged.tsv \
        -m $WORKDING_DIR/models/final_train.model.pkl \
        -p predictions/${prefix}

    $WORKDING_DIR/src/reorder.py predictions/${prefix}.prediction.txt > predictions/${prefix}.output.txt
done 


for feat in "RBP" "mmsplice" "TFBS" "eQTL" "histone" "cons-ext" "PSSM" "disorder" "AF" "chromatin-accessibility";
do
    option=drop${feat}
    echo "- "$option
    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/cv.py \
        -c $WORKDING_DIR/train/model_configs/bayes_opt_protein-snp_rfecv.json \
        -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
        -p cv_protein_snp.${option} \
        -f features_config/feature_selection.protein-snp.yml \
        --ignored $feat &> log/cv_protein_snp.${option}.log &


    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/cv.py \
        -c $WORKDING_DIR/train/model_configs/bayes_opt_protein-other_rfecv.json \
        -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
        -p cv_protein_other.${option} \
        -f features_config/feature_selection.protein-other.yml \
        --ignored $feat &> log/cv_protein_other.${option}.log &


    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/cv.py \
        -c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-snp_rfecv.json \
        -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
        -p cv_nonpro_snp.${option} \
        -f features_config/feature_selection.nonpro-snp.yml \
        --ignored $feat &> log/cv_nonpro_snp.${option}.log &


    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/cv.py \
        -c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-other_rfecv.json \
        -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
        -p cv_nonpro_other.${option} \
        -f features_config/feature_selection.nonpro-other.yml \
        --ignored $feat &> log/cv_nonpro_other.${option}.log &

    wait

    $WORKDING_DIR/src/merge_4_model.py \
        --ps predictions/cv_protein_snp.${option}.prediction.txt \
        --po predictions/cv_protein_other.${option}.prediction.txt \
        --ns predictions/cv_nonpro_snp.${option}.prediction.txt \
        --no predictions/cv_nonpro_other.${option}.prediction.txt  > predictions/train.${option}.merged.tsv
    wait

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/linearregression.py \
        -d predictions/train.${option}.merged.tsv \
        -p final_train.${option} \
        -f $WORKDING_DIR/train/feature_configs/four_group_features.yml > log/final.${option}.log
    for type in "noncoding" "coding" "missense" "synonymous";
    do 
        echo " - "${type}
        prefix=${pre}_${type}
        $WORKDING_DIR/src/split_data.py $prefix.tsv -p $prefix

        $WORKDING_DIR/src/predict.py \
            -i ${prefix}.protein_snp.tsv \
            -m $WORKDING_DIR/models/cv_protein_snp.${option}.model.pkl \
            -p predictions/${prefix}.${option}.protein_snp &> log/${prefix}.${option}.protein_snp.log &

        $WORKDING_DIR/src/predict.py \
            -i ${prefix}.protein_other.tsv \
            -m $WORKDING_DIR/models/cv_protein_other.${option}.model.pkl \
            -p predictions/${prefix}.${option}.protein_other &> log/${prefix}.${option}.protein_other.log &

        $WORKDING_DIR/src/predict.py \
            -i ${prefix}.nonpro_snp.tsv \
            -m $WORKDING_DIR/models/cv_nonpro_snp.${option}.model.pkl \
            -p predictions/${prefix}.${option}.nonpro_snp &> log/${prefix}.${option}.nonpro_snp.log &

        $WORKDING_DIR/src/predict.py \
            -i ${prefix}.nonpro_other.tsv \
            -m $WORKDING_DIR/models/cv_nonpro_other.${option}.model.pkl \
            -p predictions/${prefix}.${option}.nonpro_other &> log/${prefix}.${option}.nonpro_other.log &
        wait

        $WORKDING_DIR/src/merge_4_model.py \
            --ps predictions/${prefix}.${option}.protein_snp.prediction.txt \
            --po predictions/${prefix}.${option}.protein_other.prediction.txt \
            --ns predictions/${prefix}.${option}.nonpro_snp.prediction.txt \
            --no predictions/${prefix}.${option}.nonpro_other.prediction.txt > predictions/test_${type}.${option}.merged.tsv

        $WORKDING_DIR/src/predict_lr.py \
            -i predictions/test_${type}.${option}.merged.tsv \
            -m $WORKDING_DIR/models/final_train.${option}.model.pkl \
            -p predictions/${prefix}.${option}

        $WORKDING_DIR/src/reorder.py predictions/${prefix}.${option}.prediction.txt > predictions/${prefix}.${option}.output.txt
    done
done