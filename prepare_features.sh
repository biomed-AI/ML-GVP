#!/bin/bash

if [ $# -lt 3 ]; then
    echo "usage: $0 vcf buildver outdir"
    exit 1
fi
## input
input_vcf="$1"
export buildver="$2"
# if [ $buildver != "hg19"  ] && [ $buildver != "hg38" ]; then
if [ $buildver != "hg19"  ]; then
    echo - "Unknown build version: $buildver (hg19 required)"
    exit 1
fi
export outdir="$3"

SRC="`dirname $0`/src"

bn=`basename $input_vcf`
bn=${bn%.gz}
export bn=${bn%.vcf}

export WORKING_DIR=$(dirname `realpath $0`)
source ./config.sh

test -d $outdir || mkdir -p $outdir
feature_dir="${outdir}/features"
test -d ${feature_dir} || mkdir -p $feature_dir
temp_dir="${outdir}/temp"
test -d ${temp_dir} || mkdir -p $temp_dir




## configuration of variables
if [ $buildver = "hg19" ]; then
    export tss_bed="$WORKING_DIR/data/gencode.v38lift37.tss.bed.gz"
    export ref_fasta="$hg19"
    export gtf="$hg19_gtf"
elif [ $buildver = "hg38" ]; then
    export tss_bed=""
    export ref_fasta="$hg38"
    export gtf="$hg38_gtf"
fi

## ======================= main ============================= 
## 1. format
vcf=${outdir}/${bn}.std.vcf
$SRC/format_vcf.py $input_vcf > $vcf
$SRC/vcf_add_header.py $vcf -r $ref_fasta --overwrite
avinput=${vcf%.vcf}.avinput


## 2. ANNOVAR annotation
echo "- run annovar annotation"
${ANNOVAR_BIN}/convert2annovar.pl -format vcf4 $vcf > $avinput 2> ${avinput}.log
${ANNOVAR_BIN}/annotate_variation.pl \
    --exonicsplicing \
    --transcript_function \
    -geneanno \
    -buildver ${buildver} \
    ${avinput} \
    -dbtype ensGene \
    $HOME/db/annovar/humandb 2>> ${avinput}.log

$SRC/clean_exonic.py \
    -exonic ${avinput}.exonic_variant_function \
    -m $WORKING_DIR/data/transcript_to_uniprot.json > ${avinput}.exonic_variant_function.clean 2> ${avinput}.exonic_variant_function.clean.log || exit 1

if [ -e ${avinput}.coding_change.fa.gz ]; then
    echo "- skip coding_change"
else
    echo "- run coding_change.pl"
    ${ANNOVAR_BIN}/coding_change.pl \
        --includesnp \
        --onlyAltering \
        ${avinput}.exonic_variant_function.clean \
        ${ANNOVAR_DB}/${buildver}_ensGene.txt \
        ${ANNOVAR_DB}/${buildver}_ensGeneMrna.fa 2> ${avinput}.coding_change.log | gzip - > ${avinput}.coding_change.fa.gz
fi

$SRC/mark_coding_change.py \
    ${avinput}.coding_change.fa.gz \
    -vcf $vcf 2> ${avinput}.coding_change.marked.fa.gz.log | gzip - > ${avinput}.coding_change.marked.fa.gz

${SRC}/get_startchange.py ${avinput}.coding_change.marked.fa.gz > ${avinput}.coding_change.marked.start_change.json


## 2. prepare features
if [ $buildver == "hg19" ]; then
    hg19_vcf=$vcf
    hg38_vcf=${vcf%.vcf}.hg19To38.vcf
    CrossMap.py vcf \
        --no-comp-alleles \
        ${WORKING_DIR}/data/liftover/hg19ToHg38.over.chain.gz \
        ${vcf} $hg38 ${hg38_vcf} &> ${hg38_vcf}.log && \
    $SRC/vcf_add_header.py $hg38_vcf -r $hg38 --overwrite
else
    hg38_vcf=$vcf
    hg19_vcf=${vcf%.vcf}.hg38ToHg19.vcf
    CrossMap.py vcf \
        --no-comp-alleles \
        ${WORKING_DIR}/data/liftover/hg38ToHg19.over.chain.gz \
        ${vcf} $hg19 ${hg19_vcf} &> ${hg19_vcf}.log && \
    $SRC/vcf_add_header.py $hg19_vcf -r $hg19 --overwrite
fi


$SRC/mark_vcf.py \
  -e ${avinput}.exonic_variant_function.clean \
  -f ${avinput}.variant_function \
  -m ${WORKING_DIR}/data/transcript_to_uniprot.json \
  -tx $tss_bed $vcf > ${avinput}.marked.vcf

$SRC/vcf_add_header.py ${avinput}.marked.vcf -r $ref_fasta --overwrite


## protein features

$SRC/prepare_protein_change.py ${avinput}.marked.vcf -m ${avinput}.coding_change.marked.fa.gz > ${avinput}.features-prot-sim-change.tsv

$SRC/get_protein_subseq.py ${avinput}.coding_change.marked.fa.gz > ${avinput}.coding_change.marked.substr.fa

# run blastp
out_blastp=${avinput}.coding_change.marked.substr.uniprot.xml.gz
if [ -e ${out_blastp} ]; then
    echo "- skip $out_blastp"
else
    echo -e "run blastp => $out_blastp..."
    test -e "${uniprot_fasta}.phr" || makeblastdb -in $uniprot_fasta -dbtype prot
    $blastp \
        -outfmt 5 \
        -num_threads $MAX_THREADS \
        -evalue 0.1 \
        -num_alignments 10 \
        -query ${avinput}.coding_change.marked.substr.fa \
        -db $uniprot_fasta 2> ${avinput}.coding_change.marked.substr.uniprot.xml.gz.log | gzip - > ${avinput}.coding_change.marked.substr.uniprot.xml.gz
fi

$SRC/fetch_pssm_features.py \
    --input ${avinput}.coding_change.marked.substr.fa \
    -db $PSSM_SPOT_PATH \
    --blast-map ${avinput}.coding_change.marked.substr.uniprot.xml.gz 2> ${avinput}.pssm-disorder.pkl.gz.log && \
$SRC/prepare_protein_feats.py $vcf -p ${avinput}.coding_change.marked.substr.fa.pssm-disorder.pkl.gz 2> ${avinput}.features-pssm-disorder.tsv.log > ${avinput}.features-pssm-disorder.tsv


## hg38-only annotation
chroms=`grep -v "#" $hg38_vcf | cut -f 1 | sort | uniq`
## allele frequency
test -d ${outdir}/gnomad_af || mkdir -p ${outdir}/gnomad_af
for c in $chroms; do
    bedtools intersect -wo -a ${gnomAD_AF_DB}/gnomad.genomes.v3.1.sites.$c.AF.vcf.gz -b $hg38_vcf > ${outdir}/gnomad_af/${bn}.AF-overlap.chr${c}.hg38.tsv & done
wait
$SRC/prepare_gnomAD_AF.py ${avinput}.marked.vcf --af ${outdir}/gnomad_af/${bn}.AF-overlap.chr*.hg38.tsv > ${avinput}.features-AF.tsv


##  CRE (hg38-only annotation)
$SRC/vcf2bed.py ${hg38_vcf} > ${avinput}.hg38.bed
bedtools window -w 1000 -a ${avinput}.hg38.bed -b ${WORKING_DIR}/data/encode_cre/GRCh38-cCREs.bed | gzip - > ${avinput}.CRE.tsv.gz
$SRC/prepare_cre.py ${avinput}.marked.vcf -isec ${avinput}.CRE.tsv.gz > ${avinput}.features-CRE.tsv 2> ${avinput}.features-CRE.tsv.log

test -d ${outdir}/RBP || mkdir -p ${outdir}/RBP
for i in `seq $MAX_THREADS`; do
    # for bed in /bigdat1/pub/ENCODE/eCLIP/HepG2_AGGF1-human.bed.gz
    for bed in ${ENCODE_eCLIP_DB}/*-human.bed.gz; do
        rbp=`basename $bed -human.bed.gz`
        test -e ${bn}_${rbp}.tsv.gz && continue
        bedtools intersect -wo -a ${avinput}.hg38.bed -b $bed | gzip - > ${outdir}/RBP/${bn}.RBP_${rbp}.tsv.gz
    done &
done
wait
$SRC/prepare_eCLIP.py ${avinput}.marked.vcf --rbp ${outdir}/RBP/${bn}.RBP_*.tsv.gz > ${avinput}.features-eCLIP.tsv



## mmsplice (hg19/hg38, dependent on reference genome and annotation)
if [ -e ${avinput}.mmsplice.tsv ]; then
    echo -e "- skip mmsplice"
else
    echo -e "- run mmsplice => ${avinput}.mmsplice.tsv"
    $mmsplice_python $SRC/mmsplice_predict.py $vcf --gtf $gtf --fasta $ref_fasta -o ${avinput}.mmsplice.tsv  -bs  256 &> ${avinput}.mmsplice.tsv.log
fi
$SRC/prepare_mmsplice.py $vcf -m ${avinput}.mmsplice.tsv 2> ${avinput}.features-mmsplice.tsv.log > ${avinput}.features-mmsplice.tsv || exit "ERROR in mmsplice module!"


## 
## TFBS (hg19-only annotation)
$SRC/vcf2bed.py $vcf > ${avinput}.bed
test -d ${outdir}/chipseq || mkdir -p ${outdir}/chipseq
for bed in ${ENCODE_CHIPSeq_DB}/ENC*.hg38ToHg19.bed; do
    assay=`basename $bed .hg38ToHg19.bed`
    bedtools intersect -wo -a ${avinput}.bed -b $bed | gzip - > ${outdir}/chipseq/${bn}.${assay}-chipseq.tsv.gz &
done
wait
$SRC/prepare_chipseq.py ${avinput}.marked.vcf -isec ${outdir}/chipseq/${bn}.*-chipseq.tsv.gz > ${avinput}.features-chipseq.tsv


## stranded TFBS (from gimme, hg19-only annotation)
test -d ${outdir}/gimme || mkdir -p ${outdir}/gimme
for bed in ${WORKING_DIR}/data/encode_chipseq_gimme/*gimme.gz; do
    assay=`basename $bed .hg38ToHg19.merged.gimme.gz`
    bedtools intersect -wo -a $avinput.bed -b $bed | gzip - > ${outdir}/gimme/${bn}.${assay}-gimme.tsv.gz &
done
wait
$SRC/prepare_gimme.py ${avinput}.marked.vcf -isec ${outdir}/gimme/${bn}.*-gimme.tsv.gz > ${avinput}.features-gimme.tsv

## eQTL (hg19-only annotation)
sed 's/^/chr/;s/^chrchr/chr/' ${avinput}.bed | bedtools window -l 2000 -r 2000 -b - -a ./data/GTEx.eQTL.v8.hg38to19.merged.bed.gz >${avinput}.eQTL-2k.tsv 
$SRC/prepare_eQTL.py ${vcf} -e ${avinput}.eQTL-2k.tsv > ${avinput}.features-eQTL.tsv


## histone (hg19-only annotation)
test -d ${outdir}/roadmap_histone || mkdir -p  ${outdir}/roadmap_histone
for i in `seq $MAX_THREADS`; do
    for eid in `seq 63 113`; do
        eid=`python -c "print('E%03d' % $eid)"`

        test -e ${outdir}/roadmap_histone/${bn}.histone-${eid}-ext.tsv.gz && continue
        n=`ls ${ROADMAP_DB}/${eid}*.pval.signal.bigwig 2> /dev/null | wc -l `
        if [ $n -eq 0 ]; then
            continue
        fi
        $SRC/prepare_roadmap.py \
            ${avinput}.marked.vcf \
            --bed ${avinput}.bed \
            -bw ${ROADMAP_DB}/${eid}*.bigwig  2> ${outdir}/roadmap_histone/${bn}.histone-${eid}-ext.tsv.gz.log | gzip - > ${outdir}/roadmap_histone/${bn}.histone-${eid}-ext.tsv.gz
    done & sleep 2
done
wait
$SRC/prepare_roadmap_merge.py ${avinput}.marked.vcf -a ${outdir}/roadmap_histone/${bn}.histone-*-ext.tsv.gz > ${avinput}.features-roadmap.tsv 2> ${avinput}.features-roadmap.tsv.log


## conservation 
$SRC/prepare_cons.py \
    ${avinput}.marked.vcf \
    --bed ${avinput}.bed \
    --name phastCons100way phyloP100way phastCons46way.primates phyloP46way.primates phastCons46way.vertebrate phyloP46way.vertebrate \
    --bw \
        $hg19_phastCons_100way $hg19_phyloP_100way \
        $hg19_phastCons_46way_primates $hg19_phyloP_46way_primates \
        $hg19_phastCons_46way_vertebrate $hg19_phyloP_46way_vertebrate > ${avinput}.features-cons-ext.tsv 2> ${avinput}.features-cons-ext.tsv.log

$SRC/prepare_basic.py $avinput.marked.vcf --start ${avinput}.coding_change.marked.start_change.json  --annovar ${avinput}.exonic_variant_function.clean > ${avinput}.features-basic.tsv 2> ${avinput}.features-basic.tsv.log

$SRC/prepare_all_features.py \
    ${avinput}.features-basic.tsv \
    ${avinput}.features-CRE.tsv \
    ${avinput}.features-AF.tsv \
    ${avinput}.features-eQTL.tsv \
    ${avinput}.features-cons-ext.tsv \
    ${avinput}.features-chipseq.tsv \
    ${avinput}.features-gimme.tsv \
    ${avinput}.features-eCLIP.tsv \
    ${avinput}.features-roadmap.tsv \
    ${avinput}.features-mmsplice.tsv \
    ${avinput}.features-prot-sim-change.tsv \
    ${avinput}.features-pssm-disorder.tsv \
    -o ${avinput}.features-ALL.tsv &> ${avinput}.features-ALL.tsv.log

