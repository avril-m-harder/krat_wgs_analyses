#!/bin/bash

#SBATCH --job-name=vcf_subset
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00

module load bcftools/1.15

cd /home/amh0254/krat_roh_analyses/03_gatk/

# bcftools index krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz

### Just for samples related in some way
# bcftools view \
# --regions-file /home/amh0254/krat_roh_analyses/sample_lists/rohverlap_coordinates.txt \
# --output-type z \
# --output krat_plink_rohverlap.vcf.gz \
# krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz

# bcftools query -f '%CHROM\t%POS\t[\t%GT]\n' --output krat_plink_rohverlap_GTs.txt \
# krat_plink_rohverlap.vcf.gz


### For all pairwise sample comparisons
bcftools view \
--regions-file /home/amh0254/krat_roh_analyses/sample_lists/rohverlap_coordinates_ALL_SAMPLES.txt \
--output-type z \
--output krat_plink_rohverlap_ALL_SAMPLES.vcf.gz \
krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz

bcftools query -f '%CHROM\t%POS\t[\t%GT]\n' --output krat_plink_rohverlap_GTs_ALL_SAMPLES.txt \
krat_plink_rohverlap_ALL_SAMPLES.vcf.gz