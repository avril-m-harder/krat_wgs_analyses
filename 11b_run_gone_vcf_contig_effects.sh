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

cp /home/amh0254/krat_roh_analyses/sample_lists/go_test_sample_list.txt .
cp /home/amh0254/krat_roh_analyses/sample_lists/gone_contig_sampling/* .

# bcftools index krat_gone_filt_integerchrnames.vcf.gz


## --------------------------------
## Loop over contig sets
for k in 25 50 100
do
	mkdir /home/amh0254/krat_roh_analyses/11_gone/100_replicates_n12/${k}_contigs/
	mkdir ./${k}_contigs/
	
	bcftools view \
	--samples-file go_test_sample_list.txt \
	--output-type z \
	--output ./${k}_contigs/subset_gone_filt.recode.vcf.gz \
	--regions-file ${k}_contigs.txt \
	krat_gone_filt_integerchrnames.vcf.gz
	
	cd ./${k}_contigs/

	vcftools --gzvcf subset_gone_filt.recode.vcf.gz --plink --out gone_subset_plink
	
	sed -i "s/\t0/\t-9/g" gone_subset_plink.ped

	## --------------------------------
	## Run the program 100X 
	cp -r /home/amh0254/programs/Linux/* .

	cd ./PROGRAMMES/
	chmod +x ./*
	cd ../
	
	sed -i "s/maxNCHROM=-99/maxNCHROM=${k}/g" INPUT_PARAMETERS_FILE

	for c in {1..100}
	do
		bash script_GONE.sh gone_subset_plink
		cp Output_Ne_gone_subset_plink \
		/home/amh0254/krat_roh_analyses/11_gone/100_replicates_n12/${k}_contigs/Output_Ne_gone_subset_plink_${k}_${c}
	done
	
	cd ../
	
done

for k in 10
do
	mkdir /home/amh0254/krat_roh_analyses/11_gone/100_replicates_n12/${k}_contigs/
	mkdir ./${k}_contigs/
	
	bcftools view \
	--samples-file go_test_sample_list.txt \
	--output-type z \
	--output ./${k}_contigs/subset_gone_filt.recode.vcf.gz \
	--regions-file ${k}_contigs.txt \
	krat_gone_filt_integerchrnames.vcf.gz
	
	cd ./${k}_contigs/
	
	vcftools --gzvcf subset_gone_filt.recode.vcf.gz --plink --out gone_subset_plink
	
	sed -i "s/\t0/\t-9/g" gone_subset_plink.ped


	## --------------------------------
	## Run the program 100X 
	cp -r /home/amh0254/programs/Linux/* .

	cd ./PROGRAMMES/
	chmod +x ./*
	cd ../
	
	sed -i "s/maxNCHROM=-99/maxNCHROM=${k}/g" INPUT_PARAMETERS_FILE

	for c in {70..100}
	do
		bash script_GONE.sh gone_subset_plink
		cp Output_Ne_gone_subset_plink \
		/home/amh0254/krat_roh_analyses/11_gone/100_replicates_n12/${k}_contigs/Output_Ne_gone_subset_plink_${k}_${c}
	done
	
	cd ../
	
done

















