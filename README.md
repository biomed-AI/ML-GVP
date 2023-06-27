# ML-GVE

**ML-GVE** is an XGBoost-based model for predicting pathogenicity of variants in exons and introns. It was originally designed for the [Sherloc challenge](https://genomeinterpretation.org/cagi6-invitae.html) in the 6th Critical Assessment of Genome Interpretation. The model is one of the two best-performing models in the Sherloc challenge, ranking the first for overall and indel prediction. This repository extended the original implementation ([chenkenbio/CAGI6-Sherloc-model](https://github.com/chenkenbio/CAGI6-Sherloc-model)) and provided scripts for additional analysis and experiments.

## System Requirements and Installation
ML-GVE was trained on Linux. The key system requirements are Python (>3.8) and XGBoost (xgboost-gpu>1.6.2). The databases used for collecting features were list in [`./config.sh`]

## Bash Scripts
The configuration about the path of databases(including RoadMap, gnomAD and ENCODE, etc) and other bioinformatic tools (such as MMSplice, annovar, etc) were list in the script [`./config.sh`](https://github.com/biomed-AI/ML-GVE/blob/master/config.sh).
The script [`./prepare_features.sh`](https://github.com/biomed-AI/ML-GVE/blob/master/prepare_features.sh) can be executed to collect the features for the mutations given.
The script [`./run_predict.sh`](https://github.com/biomed-AI/ML-GVE/blob/master/run_predict.sh) can be used for predicing the pathogenicity of mutations after preparing features.

## Repository Structure and Usage
[`./train`](https://github.com/biomed-AI/ML-GVE/blob/master/train/) included all codes used for training an ML-GVE model.\
[`./src`](https://github.com/biomed-AI/ML-GVE/blob/master/src/) included all other scripts except those for model training.\
[`./caseStudy`](https://github.com/biomed-AI/ML-GVE/blob/master/caseStudy/) included the mutations from *denovo-db* dataset which was collected for case study in the paper.\
[`./cagi6-sherloc`](https://github.com/biomed-AI/ML-GVE/blob/master/cagi6-sherloc/) inlcuded all training and test data which were provided by CAGI6 competition.\
[`./ablation`](https://github.com/biomed-AI/ML-GVE/blob/master/ablation/) included those scripts used for ablation study.
