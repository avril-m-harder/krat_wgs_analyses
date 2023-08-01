#!/bin/bash

#SBATCH --job-name=ms_snps__group_
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mem=40000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set project name
PROJ=ms_kratroh_03_snp_calling

# mkdir /scratch/${USER}/${PROJ}/
chmod 700 /scratch/${USER}/${PROJ}/_group_
cd /scratch/${USER}/${PROJ}/_group_


