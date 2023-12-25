#!/bin/bash


echo "usage: $0 input"
# if [ $# -lt 1 ]; then
#     exit 1
# fi
# if [ "$2" == "all" ]
# then
#     option=""
# elif [ "$2" == "rfecv" ]
# then
#     option="_$2"
# else
#     echo "mode: null or rfecv"
#     exit 1
# fi

export WORKING_DIR=$(dirname `realpath $0`)
source ./config.sh

input="$1"
prefix=`basename $input`
prefix=${prefix%.gz}
prefix=${prefix%.tsv}

$WORKING_DIR/src/split_data.py $input -p $prefix


# $WORKING_DIR/src/predict.py -i ${prefix}.
$WORKING_DIR/src/predict.py \
    -i ${WORKING_DIR}/${prefix}.protein_snp.tsv \
    -m ${WORKING_DIR}/models/cv_protein_snp_rfecv.model.pkl \
    -p ${prefix}.protein_snp &> $WORKING_DIR/log/${prefix}.protein_snp.log

$WORKING_DIR/src/predict.py \
    -i ${WORKING_DIR}/${prefix}.protein_other.tsv \
    -m ${WORKING_DIR}/models/cv_protein_other_rfecv.model.pkl \
    -p ${prefix}.protein_other &> $WORKING_DIR/log/${prefix}.protein_other.log

$WORKING_DIR/src/predict.py \
    -i ${WORKING_DIR}/${prefix}.nonpro_snp.tsv \
    -m ${WORKING_DIR}/models/cv_nonpro_snp_rfecv.model.pkl \
    -p ${prefix}.nonpro_snp &> $WORKING_DIR/log/${prefix}.nonpro_snp.log

$WORKING_DIR/src/predict.py \
    -i ${WORKING_DIR}/${prefix}.nonpro_other.tsv \
    -m ${WORKING_DIR}/models/cv_nonpro_other_rfecv.model.pkl \
    -p ${prefix}.nonpro_other &> $WORKING_DIR/log/${prefix}.nonpro_other.log

$WORKING_DIR/src/merge_4_model.py \
    --ps ${WORKING_DIR}/${prefix}.protein_snp.prediction.txt \
    --po ${WORKING_DIR}/${prefix}.protein_other.prediction.txt \
    --ns ${WORKING_DIR}/${prefix}.nonpro_snp.prediction.txt \
    --no ${WORKING_DIR}/${prefix}.nonpro_other.prediction.txt > ${WORKING_DIR}/${prefix}.combined.tsv


$WORKING_DIR/src/predict_lr.py \
    -i ${WORKING_DIR}/${prefix}.combined.tsv \
    -m $WORKING_DIR/models/final_rfecv.model.pkl \
    -p ${prefix}.final_rfecv

$WORKING_DIR/src/reorder.py ${WORKING_DIR}/${prefix}.final_rfecv.prediction.txt > $WORKING_DIR/${prefix}_rfecv.output.txt
gzip -f $WORKING_DIR/${prefix}_rfecv.output.txt

$WORKING_DIR/train/ensemble.py \
    -i $WORKING_DIR/${prefix}_rfecv.output.txt.gz \
    -m /home/dingml/documents/CAGI6-Sherloc-model-main/models/final_ensemble.v1.model.pkl \
    -r  /home/dingml/documents/CAGI6-Sherloc-tools/data/refer_tools \
    -o $WORKING_DIR/${prefix}_rfecv.meta.output.txt
