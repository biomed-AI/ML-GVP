#!/bin/bash


echo "usage: $0 input outdir"
if [ $# -lt 1 ]; then
    exit 1
fi
export option="_rfecv"
export outdir=$3

WORKING_DIR=`dirname $0`
source ./config.sh
test -d $outdir || mkdir -p $outdir

input="$1"
prefix=`basename $input`
prefix=${prefix%.gz}
prefix=${prefix%.tsv}

$WORKING_DIR/src/split_data.py $input -p $prefix


# $WORKING_DIR/src/predict.py -i ${prefix}.
$WORKING_DIR/src/predict.py \
    -i ${prefix}.protein_snp.tsv \
    -m ${WORKING_DIR}/models/cv_protein_snp$option.model.pkl \
    -p ${prefix}.protein_snp &> $WORKING_DIR/$outdir/${prefix}.protein_snp.log

$WORKING_DIR/src/predict.py \
    -i ${prefix}.protein_other.tsv \
    -m ${WORKING_DIR}/models/cv_protein_other$option.model.pkl \
    -p ${prefix}.protein_other &> $WORKING_DIR/$outdir/${prefix}.protein_other.log

$WORKING_DIR/src/predict.py \
    -i ${prefix}.nonpro_snp.tsv \
    -m ${WORKING_DIR}/models/cv_nonpro_snp$option.model.pkl \
    -p ${prefix}.nonpro_snp &> $WORKING_DIR/$outdir/${prefix}.nonpro_snp.log

$WORKING_DIR/src/predict.py \
    -i ${prefix}.nonpro_other.tsv \
    -m ${WORKING_DIR}/models/cv_nonpro_other$option.model.pkl \
    -p ${prefix}.nonpro_other &> $WORKING_DIR/$outdir/${prefix}.nonpro_other.log

$WORKING_DIR/src/merge_4_model.py \
    --ps $WORKING_DIR/${prefix}.protein_snp.prediction.txt \
    --po $WORKING_DIR/${prefix}.protein_other.prediction.txt \
    --ns $WORKING_DIR/${prefix}.nonpro_snp.prediction.txt \
    --no $WORKING_DIR/${prefix}.nonpro_other.prediction.txt > $WORKING_DIR/${prefix}.combined.tsv


$WORKING_DIR/src/predict_lr.py \
    -i $WORKING_DIR/${prefix}.combined.tsv \
    -m $WORKING_DIR/models/final${option}.model.pkl \
    -p ${prefix}.final${option}

echo "The ML-GVP scores were provided in $WORKING_DIR/$outdir/${prefix}${option}.output.txt"
$WORKING_DIR/src/reorder.py $WORKING_DIR/${prefix}.final${option}.prediction.txt > $WORKING_DIR/$outdir/${prefix}${option}.output.txt

cat $WORKING_DIR/$outdir/${prefix}${option}.output.txt | grep -v "#" | sed 's/_/\t/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g' | sed 's/^chr//g' | sort -k1,1n -k2,2n -k3,3 -k4,4 | sed 's/^/chr/g' | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g' | awk '{printf("%s_%s_%s_%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6)}' > $WORKING_DIR/$outdir/${prefix}${option}.output.txt.mlgvp.sorted


## Ensemble
$WORKING_DIR/src/ensemble/prepare_ensemble_scores.py --ref_dir $ensemble_dir -i $WORKING_DIR/$outdir/${prefix}${option}.output.txt

# CADD
cat $WORKING_DIR/$outdir/${prefix}${option}.output.txt.CADD | sed 's/_/\t/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g' | sed 's/^chr//g' | sort -k1,1n -k2,2n -k3,3 -k4,4 | sed 's/^/chr/g' | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g' | awk '{printf("%s_%s_%s_%s\t%s\n", $1,$2,$3,$4,$5)}' > $WORKING_DIR/$outdir/${prefix}${option}.output.txt.CADD.sorted
# MVP# 
cat $WORKING_DIR/$outdir/${prefix}${option}.output.txt.MVP | sed 's/_/\t/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g' | sed 's/^chr//g' | sort -k1,1n -k2,2n -k3,3 -k4,4 | sed 's/^/chr/g' | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g' | awk '{printf("%s_%s_%s_%s\t%s\n", $1,$2,$3,$4,$5)}' > $WORKING_DIR/$outdir/${prefix}${option}.output.txt.MVP.sorted
# M-CAP
cat $WORKING_DIR/$outdir/${prefix}${option}.output.txt.MCAP | sed 's/_/\t/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g' | sed 's/^chr//g' | sort -k1,1n -k2,2n -k3,3 -k4,4 | sed 's/^/chr/g' | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g' | awk '{printf("%s_%s_%s_%s\t%s\n", $1,$2,$3,$4,$5)}' > $WORKING_DIR/$outdir/${prefix}${option}.output.txt.MCAP.sorted
# PrimateAI
cat $WORKING_DIR/$outdir/${prefix}${option}.output.txt.PrimateAI | sed 's/_/\t/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g' | sed 's/^chr//g' | sort -k1,1n -k2,2n -k3,3 -k4,4 | sed 's/^/chr/g' | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g' | awk '{printf("%s_%s_%s_%s\t%s\n", $1,$2,$3,$4,$5)}' > $WORKING_DIR/$outdir/${prefix}${option}.output.txt.PrimateAI.sorted
# AlphaMissense
cat $WORKING_DIR/$outdir/${prefix}${option}.output.txt.AlphaMissense | sed 's/_/\t/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g' | sed 's/^chr//g' | sort -k1,1n -k2,2n -k3,3 -k4,4 | sed 's/^/chr/g' | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g' | awk '{printf("%s_%s_%s_%s\t%s\n", $1,$2,$3,$4,$5)}' > $WORKING_DIR/$outdir/${prefix}${option}.output.txt.AlphaMissense.sorted


paste $WORKING_DIR/$outdir/${prefix}${option}.output.txt.mlgvp.sorted $WORKING_DIR/$outdir/${prefix}${option}.output.txt.CADD.sorted $WORKING_DIR/$outdir/${prefix}${option}.output.txt.MVP.sorted $WORKING_DIR/$outdir/${prefix}${option}.output.txt.MCAP.sorted $WORKING_DIR/$outdir/${prefix}${option}.output.txt.PrimateAI.sorted $WORKING_DIR/$outdir/${prefix}${option}.output.txt.AlphaMissense.sorted | cut -f 1-3,5,7,9,11,13| awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$4,$5,$6,$7,$8,$3)}' > $WORKING_DIR/$outdir/${prefix}${option}.combined.output.txt

$WORKING_DIR/src/ensemble/txt_addheader.py -i $WORKING_DIR/$outdir/${prefix}${option}.combined.output.txt -o $WORKING_DIR/$outdir/${prefix}${option}.combined.output.txt.extend

echo "The ML-GVP+ scores were provided in $WORKING_DIR/$outdir/${prefix}${option}.meta.output.txt"
$WORKING_DIR/src/ensemble/metapredict.py \
    -i $WORKING_DIR/$outdir/${prefix}${option}.combined.output.txt.extend \
    -o $WORKING_DIR/$outdir/${prefix}${option}.meta.output.txt \
    -m $WORKING_DIR/models/final_ensemble.v1.model.pkl
