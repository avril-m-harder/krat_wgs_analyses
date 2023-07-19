#!/bin/bash

#SBATCH --job-name=gone_resamps
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

## --------------------------------
## Copy GONE materials
cp -r /home/amh0254/programs/Linux/* .

cd ./PROGRAMMES/
chmod +x ./*
cd ../
	
	
## --------------------------------
## Iteratively copy sample list, run GONE (once per sample set)

for r in {1..100}
do

	cp /home/amh0254/krat_roh_analyses/sample_lists/gone_sample_lists/${r}.txt \
	./go_test_sample_list.txt
	
	bcftools view \
	--samples-file go_test_sample_list.txt \
	--output-type z \
	--output subset_gone_filt.recode.vcf.gz \
	krat_gone_filt_integerchrnames.vcf.gz
	
	vcftools --gzvcf subset_gone_filt.recode.vcf.gz --plink --out gone_subset_plink
	
	sed -i "s/\t0/\t-9/g" gone_subset_plink.ped

	bash script_GONE.sh gone_subset_plink
	cp Output_Ne_gone_subset_plink \
	/home/amh0254/krat_roh_analyses/11_gone/100_resamplings/Output_Ne_gone_subset_plink_${r}

done

