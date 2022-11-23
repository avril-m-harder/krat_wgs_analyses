#!/bin/bash

#SBATCH --job-name=gone
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_11_gone

## Create a directory on /scratch
mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/

## --------------------------------
## Load modules 
module load bcftools/1.11
module load vcftools/0.1.17


## --------------------------------
## Subset VCF to sample list(s) and format for PLINK.
## Final output of 11a.sh: krat_gone_filt_integerchrnames.vcf.gz (200 largest contigs
## renamed as integers 1-200)

## !!! Just trying one sample set to get things working, but can loop over lists of 
## !!! sample sets and name files accordingly to keep output from all rounds. 100 for each ## !!! of the 2 approaches should be more than plenty.

vcftools --gzvcf krat_gone_filt_integerchrnames.vcf.gz --plink --out gone_subset_plink

sed -i "s/\t0/\t-9/g" gone_subset_plink.ped


## --------------------------------
## Copy GONE materials and run the program 100X
cp -r /home/amh0254/programs/Linux/* .

cd ./PROGRAMMES/
chmod +x ./*
cd ../

for c in {1..100}
do
	bash script_GONE.sh gone_subset_plink
	cp Output_Ne_gone_subset_plink \
	/home/amh0254/krat_roh_analyses/11_gone/Output_Ne_gone_subset_plink_${c}
done


















