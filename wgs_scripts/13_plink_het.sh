#!/bin/bash

#SBATCH --job-name=plink
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 300:00:00
#SBATCH --mem=32000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_13_plink_het


## Create a directory on /scratch
mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/

# -----------------------------------------------------------------------------
# Load modules and copy input files
# -----------------------------------------------------------------------------

module load plink/1.9

cp /home/amh0254/krat_roh_analyses/06_LD_pruning/krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz .


# ----------------------------------------------------------------------
# Use PLINK to calculate heterozygosity
# ----------------------------------------------------------------------

 plink \
 	--het small-sample \
 	--allow-extra-chr \
 	--vcf krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz


