#!/bin/bash

#SBATCH --job-name=vcf_depth_check
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00


## -----------------------------------
## Load modules
module load vcftools/0.1.17
module load bcftools/1.15

#cd /home/amh0254/krat_roh_analyses/03_gatk/


## -----------------------------------
## Run vcftools to output mean depth for each individual
# vcftools --gzvcf krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz \
# --out krat_mean_indiv_depths \
# --depth

# bcftools stats -s - krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz


cd /home/amh0254/krat_roh_analyses/06_LD_pruning/

bcftools stats -s - krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz > \
krat_allcontigs_LDpruned_stats.txt


