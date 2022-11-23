#!/bin/bash

#SBATCH --job-name=_group__gtype
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mem=60000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_03_gatk

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/

## --------------------------------
## Load modules 
module load bcftools/1.11
module load vcftools/0.1.17


## --------------------------------
## Concatenate final VCF files of SNPs
# ls final_filtered_JAHHPX01000*.vcf > vcf_file_list.txt

# bcftools concat -f vcf_file_list.txt \
# -o final_allfiltcontigs_all_samps.vcf \
# --threads 6

# vcftools --vcf final_allfiltcontigs_all_samps.vcf \
# --recode --recode-INFO-all \
# --out krat_final_final_allfiltcontigs_all_samps \
# --maf 0.05 \
# --max-missing 0.75 \
# --minQ 20

gzip -c krat_final_final_allfiltcontigs_all_samps.recode.vcf > \
krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz