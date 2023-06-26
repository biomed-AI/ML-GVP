# ML-GVE
ML-GVE is an XGBoost-based model for predicting pathogenicity of variants in exons and introns. It was originally designed for the [Sherloc challenge](https://genomeinterpretation.org/cagi6-invitae.html) in the 6th Critical Assessment of Genome Interpretation. The model is one of the two best-performing models in the Sherloc challenge, ranking the first for overall and indel prediction. This repository extended the original implementation ([chenkenbio/CAGI6-Sherloc-model](https://github.com/chenkenbio/CAGI6-Sherloc-model)) and provided scripts for additional analysis and experiments.

## bash scripts
The configuration about the path of database and other tools (such as MMSplice etc) were list in the script *config.sh*.
The script *prepare_features.sh* were executed for preparing the features for the mutations given.
The script *run_predict.sh* were used for predicing the pathogenicity of mutations after preparing features.

## Folds
### train
The folder **train** included all codes used for training a ML-GVE model. 


### src
The folder **src** included all other scipts except those for model training.


### caseStudy
The folder **caseStudy** included the mutations from denovo-db dataset which was used for caseStudy in the paper.

### cagi6-sherloc
The folder **cagi6-sherloc** included all training the test data which were provided by CAGI6.


### ablation
The folder **ablation** included those bash scripts which were used for ablation study.
