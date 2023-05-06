#!/bin/bash

if [ $# -lt 1 ]; then
    echo "usage: $0 input"
    exit 1
fi

WORKDING_DIR=`dirname $0`
WORKDING_DIR=$WORKDING_DIR/../

prefix=`basename $1`
prefix=${prefix%.tsv.gz}

features=(basic RBP mmsplice TFBS eQTL histone cons-ext PSSM disorder AF chromatin-accessibility)


features_back=("${features[@]}")


$WORKDING_DIR/src/split_data.py $1 -p $prefix
for ((i=0;i<${#features[@]};i++));
do 
    suffix=${features[i]}
    unset features[i]
    # echo ${features[@]}
    echo "- "$suffix

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/cv.py \
    -c $WORKDING_DIR/train/model_configs/bayes_opt_protein-snp_rfecv.json \
    -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
    -p cv_protein_snp.${suffix} \
    -f features_config/all_features.yml \
    --ignored ${features[@]} basic &> log/cv_protein_snp.${suffix}.log &
    
    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/cv.py \
    -c $WORKDING_DIR/train/model_configs/bayes_opt_protein-other_rfecv.json \
    -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
    -p cv_protein_other.${suffix} \
    -f features_config/all_features.yml \
    --ignored ${features[@]} basic &> log/cv_protein_other.${suffix}.log &

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/cv.py \
    -c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-snp_rfecv.json \
    -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
    -p cv_nonpro_snp.${suffix} \
    -f features_config/no_protein_features.yml \
    --ignored ${features[@]} basic &> log/cv_nonpro_snp.${suffix}.log &

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/cv.py \
    -c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-other_rfecv.json \
    -d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
    -p cv_nonpro_other.${suffix} \
    -f features_config/no_protein_features.yml \
    --ignored ${features[@]} basic &> log/cv_nonpro_other.${suffix}.log &
    wait

    $WORKDING_DIR/src/merge_4_model.py \
        --ps predictions/cv_protein_snp.${suffix}.prediction.txt \
        --po predictions/cv_protein_other.${suffix}.prediction.txt \
        --ns predictions/cv_nonpro_snp.${suffix}.prediction.txt \
        --no predictions/cv_nonpro_other.${suffix}.prediction.txt  > predictions/train.${suffix}.merged.tsv

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/train/linearregression.py \
        -d predictions/train.${suffix}.merged.tsv \
        -p final_train.${suffix} \
        -f $WORKDING_DIR/train/feature_configs/four_group_features.yml > log/final.${suffix}.log
    wait

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/src/predict.py \
        -i ${WORKDING_DIR}/${prefix}.protein_snp.tsv \
        -m $WORKDING_DIR/models/cv_protein_snp.${suffix}.model.pkl \
        -p predictions/${prefix}.${suffix}.protein_snp &> log/${prefix}.${suffix}.protein_snp.log &

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/src/predict.py \
        -i ${WORKDING_DIR}/${prefix}.protein_other.tsv \
        -m $WORKDING_DIR/models/cv_protein_other.${suffix}.model.pkl \
        -p predictions/${prefix}.${suffix}.protein_other &> log/${prefix}.${suffix}.protein_other.log &

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/src/predict.py \
        -i ${WORKDING_DIR}/${prefix}.nonpro_snp.tsv \
        -m $WORKDING_DIR/models/cv_nonpro_snp.${suffix}.model.pkl \
        -p predictions/${prefix}.${suffix}.nonpro_snp &> log/${prefix}.${suffix}.nonpro_snp.log &

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/src/predict.py \
        -i ${WORKDING_DIR}/${prefix}.nonpro_other.tsv \
        -m $WORKDING_DIR/models/cv_nonpro_other.${suffix}.model.pkl \
        -p predictions/${prefix}.${suffix}.nonpro_other &> log/${prefix}.${suffix}.nonpro_other.log &
    wait

    $WORKDING_DIR/src/merge_4_model.py \
        --ps predictions/${prefix}.${suffix}.protein_snp.prediction.txt \
        --po predictions/${prefix}.${suffix}.protein_other.prediction.txt \
        --ns predictions/${prefix}.${suffix}.nonpro_snp.prediction.txt \
        --no predictions/${prefix}.${suffix}.nonpro_other.prediction.txt > predictions/test_${suffix}.merged.tsv

    export CUDA_VISIBLE_DEVICES=3
    $WORKDING_DIR/src/predict_lr.py \
        -i predictions/test_${suffix}.merged.tsv \
        -m $WORKDING_DIR/models/final_train.${suffix}.model.pkl \
        -p predictions/${prefix}.${suffix}

    $WORKDING_DIR/src/reorder.py predictions/${prefix}.${suffix}.prediction.txt > predictions/${prefix}.${suffix}.output.txt


    features=("${features_back[@]}")
done

# Others
suffix="Others"
echo "- "$suffix
export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_protein-snp_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
-p cv_protein_snp.${suffix} \
-f features_config/all_features.yml \
--ignored  RBP mmsplice  PSSM disorder  &> log/cv_protein_snp.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_protein-other_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
-p cv_protein_other.${suffix} \
-f features_config/all_features.yml \
--ignored RBP mmsplice  PSSM disorder &> log/cv_protein_other.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-snp_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
-p cv_nonpro_snp.${suffix} \
-f features_config/no_protein_features.yml \
--ignored RBP mmsplice  PSSM disorder &> log/cv_nonpro_snp.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-other_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
-p cv_nonpro_other.${suffix} \
-f features_config/no_protein_features.yml \
--ignored RBP mmsplice  PSSM disorder &> log/cv_nonpro_other.${suffix}.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/cv_protein_snp.${suffix}.prediction.txt \
    --po predictions/cv_protein_other.${suffix}.prediction.txt \
    --ns predictions/cv_nonpro_snp.${suffix}.prediction.txt \
    --no predictions/cv_nonpro_other.${suffix}.prediction.txt  > predictions/train.${suffix}.merged.tsv

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/linearregression.py \
    -d predictions/train.${suffix}.merged.tsv \
    -p final_train.${suffix} \
    -f $WORKDING_DIR/train/feature_configs/four_group_features.yml > log/final.${suffix}.log
wait

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.protein_snp.tsv \
    -m $WORKDING_DIR/models/cv_protein_snp.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.protein_snp &> log/${prefix}.${suffix}.protein_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.protein_other.tsv \
    -m $WORKDING_DIR/models/cv_protein_other.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.protein_other &> log/${prefix}.${suffix}.protein_other.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.nonpro_snp.tsv \
    -m $WORKDING_DIR/models/cv_nonpro_snp.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.nonpro_snp &> log/${prefix}.${suffix}.nonpro_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.nonpro_other.tsv \
    -m $WORKDING_DIR/models/cv_nonpro_other.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.nonpro_other &> log/${prefix}.${suffix}.nonpro_other.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/${prefix}.${suffix}.protein_snp.prediction.txt \
    --po predictions/${prefix}.${suffix}.protein_other.prediction.txt \
    --ns predictions/${prefix}.${suffix}.nonpro_snp.prediction.txt \
    --no predictions/${prefix}.${suffix}.nonpro_other.prediction.txt > predictions/test_${suffix}.merged.tsv

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict_lr.py \
    -i predictions/test_${suffix}.merged.tsv \
    -m $WORKDING_DIR/models/final_train.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}

$WORKDING_DIR/src/reorder.py predictions/${prefix}.${suffix}.prediction.txt > predictions/${prefix}.${suffix}.output.txt

# RNA
suffix="RNA"
echo "- "$suffix
export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_protein-snp_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
-p cv_protein_snp.${suffix} \
-f features_config/all_features.yml \
--ignored basic TFBS eQTL histone cons-ext PSSM disorder AF chromatin-accessibility &> log/cv_protein_snp.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_protein-other_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
-p cv_protein_other.${suffix} \
-f features_config/all_features.yml \
--ignored basic TFBS eQTL histone cons-ext PSSM disorder AF chromatin-accessibility &> log/cv_protein_other.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-snp_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
-p cv_nonpro_snp.${suffix} \
-f features_config/no_protein_features.yml \
--ignored basic TFBS eQTL histone cons-ext PSSM disorder AF chromatin-accessibility &> log/cv_nonpro_snp.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-other_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
-p cv_nonpro_other.${suffix} \
-f features_config/no_protein_features.yml \
--ignored basic TFBS eQTL histone cons-ext PSSM disorder AF chromatin-accessibility &> log/cv_nonpro_other.${suffix}.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/cv_protein_snp.${suffix}.prediction.txt \
    --po predictions/cv_protein_other.${suffix}.prediction.txt \
    --ns predictions/cv_nonpro_snp.${suffix}.prediction.txt \
    --no predictions/cv_nonpro_other.${suffix}.prediction.txt  > predictions/train.${suffix}.merged.tsv

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/linearregression.py \
    -d predictions/train.${suffix}.merged.tsv \
    -p final_train.${suffix} \
    -f $WORKDING_DIR/train/feature_configs/four_group_features.yml > log/final.${suffix}.log
wait

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.protein_snp.tsv \
    -m $WORKDING_DIR/models/cv_protein_snp.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.protein_snp &> log/${prefix}.${suffix}.protein_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.protein_other.tsv \
    -m $WORKDING_DIR/models/cv_protein_other.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.protein_other &> log/${prefix}.${suffix}.protein_other.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.nonpro_snp.tsv \
    -m $WORKDING_DIR/models/cv_nonpro_snp.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.nonpro_snp &> log/${prefix}.${suffix}.nonpro_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.nonpro_other.tsv \
    -m $WORKDING_DIR/models/cv_nonpro_other.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.nonpro_other &> log/${prefix}.${suffix}.nonpro_other.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/${prefix}.${suffix}.protein_snp.prediction.txt \
    --po predictions/${prefix}.${suffix}.protein_other.prediction.txt \
    --ns predictions/${prefix}.${suffix}.nonpro_snp.prediction.txt \
    --no predictions/${prefix}.${suffix}.nonpro_other.prediction.txt > predictions/test_${suffix}.merged.tsv

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict_lr.py \
    -i predictions/test_${suffix}.merged.tsv \
    -m $WORKDING_DIR/models/final_train.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}

$WORKDING_DIR/src/reorder.py predictions/${prefix}.${suffix}.prediction.txt > predictions/${prefix}.${suffix}.output.txt

# Protein

suffix="protein"
echo "- "$suffix
export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_protein-snp_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-snp.tsv \
-p cv_protein_snp.${suffix} \
-f features_config/all_features.yml \
--ignored basic RBP mmsplice TFBS eQTL histone cons-ext AF chromatin-accessibility &> log/cv_protein_snp.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_protein-other_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.protein-other.tsv \
-p cv_protein_other.${suffix} \
-f features_config/all_features.yml \
--ignored basic RBP mmsplice TFBS eQTL histone cons-ext AF chromatin-accessibility &> log/cv_protein_other.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-snp_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-snp.tsv \
-p cv_nonpro_snp.${suffix} \
-f features_config/no_protein_features.yml \
--ignored basic RBP mmsplice TFBS eQTL histone cons-ext AF chromatin-accessibility &> log/cv_nonpro_snp.${suffix}.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/cv.py \
-c $WORKDING_DIR/train/model_configs/bayes_opt_nonpro-other_rfecv.json \
-d $WORKDING_DIR/train/features/final_train_092921.raw.std.avinput.features-ALL.nonpro-other.tsv \
-p cv_nonpro_other.${suffix} \
-f features_config/no_protein_features.yml \
--ignored basic RBP mmsplice TFBS eQTL histone cons-ext AF chromatin-accessibility &> log/cv_nonpro_other.${suffix}.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/cv_protein_snp.${suffix}.prediction.txt \
    --po predictions/cv_protein_other.${suffix}.prediction.txt \
    --ns predictions/cv_nonpro_snp.${suffix}.prediction.txt \
    --no predictions/cv_nonpro_other.${suffix}.prediction.txt  > predictions/train.${suffix}.merged.tsv

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/linearregression.py \
    -d predictions/train.${suffix}.merged.tsv \
    -p final_train.${suffix} \
    -f $WORKDING_DIR/train/feature_configs/four_group_features.yml > log/final.${suffix}.log
wait

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.protein_snp.tsv \
    -m $WORKDING_DIR/models/cv_protein_snp.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.protein_snp &> log/${prefix}.${suffix}.protein_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.protein_other.tsv \
    -m $WORKDING_DIR/models/cv_protein_other.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.protein_other &> log/${prefix}.${suffix}.protein_other.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.nonpro_snp.tsv \
    -m $WORKDING_DIR/models/cv_nonpro_snp.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.nonpro_snp &> log/${prefix}.${suffix}.nonpro_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.nonpro_other.tsv \
    -m $WORKDING_DIR/models/cv_nonpro_other.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}.nonpro_other &> log/${prefix}.${suffix}.nonpro_other.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/${prefix}.${suffix}.protein_snp.prediction.txt \
    --po predictions/${prefix}.${suffix}.protein_other.prediction.txt \
    --ns predictions/${prefix}.${suffix}.nonpro_snp.prediction.txt \
    --no predictions/${prefix}.${suffix}.nonpro_other.prediction.txt > predictions/test_${suffix}.merged.tsv

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict_lr.py \
    -i predictions/test_${suffix}.merged.tsv \
    -m $WORKDING_DIR/models/final_train.${suffix}.model.pkl \
    -p predictions/${prefix}.${suffix}

$WORKDING_DIR/src/reorder.py predictions/${prefix}.${suffix}.prediction.txt > predictions/${prefix}.${suffix}.output.txt


# overall
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
-f features_config/feature_selection.nonpro-other.yml  &> log/cv_nonpro_other.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/cv_protein_snp.prediction.txt \
    --po predictions/cv_protein_other.prediction.txt \
    --ns predictions/cv_nonpro_snp.prediction.txt \
    --no predictions/cv_nonpro_other.prediction.txt  > predictions/train.merged.tsv

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/train/linearregression.py \
    -d predictions/train.merged.tsv \
    -p final_train \
    -f $WORKDING_DIR/train/feature_configs/four_group_features.yml > log/final.log
wait

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.protein_snp.tsv \
    -m $WORKDING_DIR/models/cv_protein_snp.model.pkl \
    -p predictions/${prefix}.protein_snp &> log/${prefix}protein_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.protein_other.tsv \
    -m $WORKDING_DIR/models/cv_protein_other.model.pkl \
    -p predictions/${prefix}.protein_other &> log/${prefix}.protein_other.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.nonpro_snp.tsv \
    -m $WORKDING_DIR/models/cv_nonpro_snp.model.pkl \
    -p predictions/${prefix}.nonpro_snp &> log/${prefix}.nonpro_snp.log &

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict.py \
    -i ${WORKDING_DIR}/${prefix}.nonpro_other.tsv \
    -m $WORKDING_DIR/models/cv_nonpro_other.model.pkl \
    -p predictions/${prefix}.nonpro_other &> log/${prefix}.nonpro_other.log &
wait

$WORKDING_DIR/src/merge_4_model.py \
    --ps predictions/${prefix}.protein_snp.prediction.txt \
    --po predictions/${prefix}.protein_other.prediction.txt \
    --ns predictions/${prefix}.nonpro_snp.prediction.txt \
    --no predictions/${prefix}.nonpro_other.prediction.txt > predictions/test.merged.tsv

export CUDA_VISIBLE_DEVICES=3
$WORKDING_DIR/src/predict_lr.py \
    -i predictions/test.merged.tsv \
    -m $WORKDING_DIR/models/final_train.model.pkl \
    -p predictions/${prefix}

$WORKDING_DIR/src/reorder.py predictions/${prefix}.prediction.txt > predictions/${prefix}.output.txt
