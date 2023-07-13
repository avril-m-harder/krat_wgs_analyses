#!/bin/bash

#SBATCH --job-name=ld_prune
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16000
#SBATCH -t 300:00:00


##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_06_LD_pruning

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
# chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/


# -----------------------------------------------------------------------------
# Load modules and copy files
# -----------------------------------------------------------------------------
module load vcftools/0.1.17
module load plink/1.9
module load bcftools/1.15

# cp /home/amh0254/krat_roh_analyses/03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz .

# cp /home/amh0254/krat_roh_analyses/03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf .


# -----------------------------------------------------------------------------
# Make chromosome map
# -----------------------------------------------------------------------------
#### the below code doesn't work (awk error)
# bcftools view -H krat_final_final_allfiltcontigs_all_samps.recode.vcf | cut -f 1 | \
# uniq > names.txt 

# awk '{print $0"\t"$0}' names.txt > krat_filtcontigs.chrom-map.txt


# -----------------------------------------------------------------------------
# Run LD-pruning in PLINK
# -----------------------------------------------------------------------------
## Produce .map file
# vcftools --plink \
# --out krat_filtcontigs \
# --chrom-map krat_filtcontigs.chrom-map.txt \
# --gzvcf krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz


## ID SNPs in strong LD -- requires 16 Gb memory (on ASC)
# plink \
# --file krat_filtcontigs \
# --out krat_filtcontigs_LDpruned \
# --allow-extra-chr \
# --indep-pairwise 50 5 0.5

## Produce a new VCF without LD-pruned SNPs
# sed 's/\:/\t/g' krat_filtcontigs_LDpruned.prune.in > LDpruned_keepsites.tsv

# vcftools \
# --gzvcf krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz \
# --out krat_final_final_allfiltcontigs_all_samps_LDpruned \
# --positions LDpruned_keepsites.tsv \
# --recode --recode-INFO-all


# -----------------------------------------------------------------------------
# Copy output back to /home/
# -----------------------------------------------------------------------------
# gzip krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf
# 
# cp krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz \
# /home/amh0254/krat_roh_analyses/06_LD_pruning/

# -----------------------------------------------------------------------------
# 11/21/22: adding in filtering for est-sfs to allow 0 missingness
# -----------------------------------------------------------------------------
# vcftools --gzvcf krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz \
# --out krat_LDpruned_zeromissingness \
# --max-missing 1 \
# --recode --recode-INFO-all

gzip krat_LDpruned_zeromissingness.recode.vcf

cp krat_LDpruned_zeromissingness.recode.vcf.gz \
/home/amh0254/krat_roh_analyses/06_LD_pruning/




