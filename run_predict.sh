#!/bin/bash


echo "usage: $0 input mode(null, rfecv)"
if [ $# -lt 1 ]; then
    exit 1
fi
if [ ! "$2" ]
then
    option=$2
elif [ "$2" == "rfecv" ]
then
    option="_$2"
else
    echo "mode: null or rfecv"
    exit 1
fi
WORKING_DIR=`dirname $0`


input="$1"
prefix=`basename $input`
prefix=${prefix%.gz}
prefix=${prefix%.tsv}

$WORKING_DIR/src/split_data.py $input -p $prefix



# $WORKING_DIR/src/predict.py -i ${prefix}.
$WORKING_DIR/src/predict.py \
    -i ${prefix}.protein_snp.tsv \
    -m ${WORKING_DIR}/models/cv_protein_snp$option.model.pkl \
    -p ${prefix}.protein_snp &> $WORKING_DIR/log/${prefix}.protein_snp.log

$WORKING_DIR/src/predict.py \
    -i ${prefix}.protein_other.tsv \
    -m ${WORKING_DIR}/models/cv_protein_other$option.model.pkl \
    -p ${prefix}.protein_other &> $WORKING_DIR/log/${prefix}.protein_other.log

$WORKING_DIR/src/predict.py \
    -i ${prefix}.nonpro_snp.tsv \
    -m ${WORKING_DIR}/models/cv_nonpro_snp$option.model.pkl \
    -p ${prefix}.nonpro_snp &> $WORKING_DIR/log/${prefix}.nonpro_snp.log

$WORKING_DIR/src/predict.py \
    -i ${prefix}.nonpro_other.tsv \
    -m ${WORKING_DIR}/models/cv_nonpro_other$option.model.pkl \
    -p ${prefix}.nonpro_other &> $WORKING_DIR/log/${prefix}.nonpro_other.log

$WORKING_DIR/src/merge_4_model.py \
    --ps $WORKING_DIR/${prefix}.protein_snp.prediction.txt \
    --po $WORKING_DIR/${prefix}.protein_other.prediction.txt \
    --ns $WORKING_DIR/${prefix}.nonpro_snp.prediction.txt \
    --no $WORKING_DIR/${prefix}.nonpro_other.prediction.txt > $WORKING_DIR/${prefix}.combined.tsv


$WORKING_DIR/src/predict_lr.py \
    -i $WORKING_DIR/${prefix}.combined.tsv \
    -m $WORKING_DIR/models/final${option}.model.pkl \
    -p ${prefix}.final${option}

$WORKING_DIR/src/reorder.py $WORKING_DIR/${prefix}.final${option}.prediction.txt > $WORKING_DIR/${prefix}${option}.output.txt