# Multi-level feature Genomic Variants Estimator (ML-GVE)

**ML-GVE** is an XGBoost-based model for predicting pathogenicity of variants in exons and introns. It was originally designed for the [Sherloc challenge](https://genomeinterpretation.org/cagi6-invitae.html) in the 6th Critical Assessment of Genome Interpretation. The model is one of the two best-performing models in the Sherloc challenge, ranking the first for overall and indel prediction. This repository extended the original implementation ([chenkenbio/CAGI6-Sherloc-model](https://github.com/chenkenbio/CAGI6-Sherloc-model)) and provided scripts for additional analysis and experiments.

## System Requirements and Installation
ML-GVE was trained on Linux. The key system requirements are Python (>3.8) and XGBoost (xgboost-gpu>1.6.2). The databases used for collecting features were list in [`./config.sh`](https://github.com/biomed-AI/ML-GVE/blob/master/config.sh).
 - [`ANNOVAR`](https://annovar.openbioinformatics.org/en/latest/)
 - [`blast: NCBI BLAST`](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
 - [`MMSplice`](https://github.com/gagneurlab/MMSplice_MTSplice)
 - [`bedtools`](https://bedtools.readthedocs.io/en/latest/)
 - [`pyBigWig`](https://github.com/deeptools/pyBigWig)
 - [`XGBoost`](https://xgboost.readthedocs.io/en/stable/)
 - [`CrossMap`](https://crossmap.sourceforge.net/#installation)
 - [`skopt`](https://scikit-optimize.github.io/)(required only during training process)
 - [`tabix`](https://pypi.org/project/pytabix/) (required to builed the .tbi file for ensemble predictors)

## Configure config.sh
Set the paths in [`./config.sh`](https://github.com/biomed-AI/ML-GVE/blob/master/config.sh):
 - `blastp`: path of the tool [`blast: NCBI BLAST`](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
 - `hg19/hg38`: GRCh37(hg19)/GRCh38(hg38) reference genome, can be visited from [`here`](https://grch37.ensembl.org/index.html) 
 - `hg19_gtf`: hg19 gene annotation in GTF format
 - `PSSM_SPOT_PATH`: path to PSSM_SPOT_DISORDER folder (decompress PSSM-HHM-SPOT_DISORDER.tgx in [`data/PSSM-HHM-SPOT_DISORDER`]
 - `ROADMAP_DB`: path to RoadMap ChIP-seq/DNase-seq processed data(download using the following commands in data/roadmap)

```
for eid in `seq 63 113`; do
    eid=`python -c "print('E%03d' % $eid)"`
    for assay in DNase H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3; do
        wget -c https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/${eid}-${assay}.pval.signal.bigwig
    done
done
```

 - `ENCODE_CHIPSeq_DB`: path to data/encode_chipseq folder (unzip data/encode_chipseq/encode_chipseq.tgz)
 - `ENCODE_eCLIP_DB`: path to RNA binding protein (RBP) binding sites data, tar -xf ENCODE_RBP.tar in data/eCLIP
 - `gtex_eQTL_loci`: path to data/GTEx.eQTL.v8.hg38to19.merged.bed.gz
 - `gnomAD_AF_DB`: path to gnomAD.v3 folder (tar -xf gnomad.v3.AF.vcf.tar in data/gnomAD.v3)
 - `uniprot_fasta`: uniprot protein fasta file (unzip data/uniprot_human.fa.gz)
 - `hg19_phastCons_100way`: [`download`](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons100way/hg19.100way.phastCons.bw)
 - `hg19_phyloP_100way`: [`download`](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/hg19.100way.phyloP100way.bw)
 - `hg19_phastCons_46way_primates`: [`download`](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates.phastCons46way.bw)
 - `hg19_phyloP_46way_primates`: [`download`](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates.phyloP46way.bw)
 - `hg19_phastCons_46way_vertebrate`: [`download`](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate.phastCons46way.bw)
 - `hg19_phyloP_46way_vertebrate`: [`download`](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate.phyloP46way.bw) \
Or download the above files in bigWig format in data/ucsc/hg19 using

```
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate.phyloP46way.bw
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates.phyloP46way.bw
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates.phastCons46way.bw
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate.phastCons46way.bw
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons100way/hg19.100way.phastCons.bw  
```
 - MMSplice: install mmsplice and set `mmsplice_python` to the python interpreter with MMSplice installed
 - ensemble_dir: the directory including CADD, MVP, M-CAP, PrimateAI, AlphaMissense (BED file and .tbi file)

## Demo
1. Prepare features:
```
./prepare_features.sh input.vcf hg19 output
```
After executing the command, the file about features can be found at `output/input.std.avinput.features-ALL.tsv`
2. Make prediction
```
./run_prediction.sh output/input.stc.avinput.features-ALL.tsv predict_dir
```
The prediction scores can be found at `input.std.avinput.features-ALL_rfecv.output.txt`. The scores of ML-GVE(meta) can be found at `input.std.avinput.features-ALL_rfecv.meta.output.txt`

## Repository Structure and Usage
[`./train`](https://github.com/biomed-AI/ML-GVE/blob/master/train/) included all codes used for training an ML-GVE model.\
[`./src`](https://github.com/biomed-AI/ML-GVE/blob/master/src/) included all other scripts except those for model training.\
[`./cagi6-sherloc`](https://github.com/biomed-AI/ML-GVE/blob/master/cagi6-sherloc/) inlcuded all training and test data which were provided by CAGI6 competition.
