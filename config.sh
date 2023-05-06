#!/bin/bash

export MAX_THREADS=16

export ANNOVAR_BIN="/home/chenken/tools/annovar"
export ANNOVAR_DB="/home/chenken/db/annovar/humandb"
export blastp="/home/chenken/opt/miniconda3/envs/blast/bin/blastp"


## reference dataset
export hg19="/home/chenken/db/gencode/GRCh37/GRCh37.primary_assembly.genome.fa"
export hg38="/home/chenken/db/gencode/GRCh38/GRCh38.primary_assembly.genome.fa"
export hg19_gtf="/bigdat1/pub/gencode/GRCh37/release_34/gencode.v34lift37.annotation.gtf"

export PSSM_SPOT_PATH="$WORKING_DIR/data/PSSM-HHM-SPOT_DISORDER"

export ROADMAP_DB="$WORKING_DIR/data/ChIP-seq+DNase-seq"

export ENCODE_CHIPSeq_DB="$WORKING_DIR/data/encode-chipseq"
export ENCODE_eCLIP_DB="$WORKING_DIR/data/eCLIP/old/for_cagi6"
export gtex_eQTL_loci="$WORKING_DIR/data/GTEx.eQTL.v8.hg38to19.merged.bed.gz"
export gnomAD_AF_DB="$WORKING_DIR/data/gnomAD.v3"

export uniprot_fasta="$WORKING_DIR/data/uniprot_human.fa"

export hg19_phastCons_100way="$WORKING_DIR/data/ucsc_cons/hg19.100way.phastCons.bw"
export hg19_phyloP_100way="$WORKING_DIR/data/ucsc/hg19/hg19.100way.phyloP100way.bw"
export hg19_phastCons_46way_primates="$WORKING_DIR/data/ucsc/hg19/primates.phastCons46way.bw"
export hg19_phyloP_46way_primates="$WORKING_DIR/data/ucsc/hg19/primates.phyloP46way.bw"
export hg19_phastCons_46way_vertebrate="$WORKING_DIR/data/ucsc/hg19/vertebrate.phastCons46way.bw"
export hg19_phyloP_46way_vertebrate="$WORKING_DIR/data/ucsc/hg19/vertebrate.phyloP46way.bw"


## setup mmsplice
# use GPU
export LD_LIBRARY_PATH="/home/chenken/opt/miniconda3/envs/mmsplice-gpu/lib:$LD_LIBRARY_PATH"
export mmsplice_python=/home/chenken/opt/miniconda3/envs/mmsplice-gpu/bin/python3
# export LD_LIBRARY_PATH="/home/dingml/anaconda3/envs/mmsplic/lib:$LD_LIBRARY_PATH"
# export mmsplice_python=/home/dingml/anaconda3/envs/mmsplice/bin/python3 
## or use CPU
#export mmsplice_python=`which python3`  # path to python with mmsplice installed
