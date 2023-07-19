#!/bin/bash

#SBATCH --job-name=vcf_copy
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

cd /scratch/avrilh/kratroh_03_gatk
cp ./*.g.vcf* /home/amh0254/krat_roh_analyses/03_gatk/
