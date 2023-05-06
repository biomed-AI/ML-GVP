#!/usr/bin/env python3
import argparse, os

from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table#, predict_save
from mmsplice.utils import max_varEff

from utils import idle_gpu
# os.environ["CUDA_VISIBLE_DEVICES"] = str(idle_gpu())
os.environ["CUDA_VISIBLE_DEVICES"] = str(idle_gpu())

import tensorflow as tf
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Output interpretation: https://github.com/gagneurlab/MMSplice_MTSplice#output")
    p.add_argument("vcf", help="VCF, should be in standard VCF format with header (A demo VCF header is in the end of the script)")
    p.add_argument("--gtf", required=True, help="GTF file: gencode.v34lift37.annotation.gtf")
    p.add_argument("--fasta", help="GRCh37.primary_assembly.genome.fa")
    p.add_argument("-bs", "--batch_size", help="batch size of MMSplice", default=32, type=int)
    p.add_argument('-o', "--output", required=True, help="Output csv file")
    return p.parse_args()

if __name__ == "__main__":
    print('\n')
    args = get_args()
    vcf = args.vcf
    gtf = args.gtf
    fasta = args.fasta
    csv = args.output
    bs = args.batch_size

    # dataloader to load variants from vcf
    dl = SplicingVCFDataloader(gtf, fasta, vcf)

    # Specify model
    model = MMSplice()

    # predict and save to csv file
    #predict_save(model, dl, csv, pathogenicity=True, splicing_efficiency=True)

    # Or predict and return as df
    predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True, batch_size=bs)

    # Summerize with maximum effect size
    predictionsMax = max_varEff(predictions)

    predictionsMax.to_csv(args.output)



#>>====================== VCF header =============================<<#

##fileformat=VCFv4.0
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr4,length=191154276>
##contig=<ID=chr5,length=180915260>
##contig=<ID=chr6,length=171115067>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chr8,length=146364022>
##contig=<ID=chr9,length=141213431>
##contig=<ID=chr10,length=135534747>
##contig=<ID=chr11,length=135006516>
##contig=<ID=chr12,length=133851895>
##contig=<ID=chr13,length=115169878>
##contig=<ID=chr14,length=107349540>
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr16,length=90354753>
##contig=<ID=chr17,length=81195210>
##contig=<ID=chr18,length=78077248>
##contig=<ID=chr19,length=59128983>
##contig=<ID=chr20,length=63025520>
##contig=<ID=chr21,length=48129895>
##contig=<ID=chr22,length=51304566>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chrY,length=59373566>
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

#>>===============================================================<<#
