#!/bin/bash

#SBATCH --job-name=get_GTs
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00

module load bcftools/1.15

cd /home/amh0254/krat_roh_analyses/03_gatk/

bcftools query -f '%CHROM\t%POS\t[\t%GT]\n' \
--output krat_LDpruned_extracted_genotypes.txt \
krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz