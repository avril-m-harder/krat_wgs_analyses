#!/bin/bash

#SBATCH --job-name=_group__gtype
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=60000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_03_gatk

## Set batch group
GROUP=_group_

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
# chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/


## --------------------------------
## Load modules 
module load picard/2.23.9
module load gatk/4.1.9.0
module load bcftools/1.11
module load samtools/1.11


## --------------------------------
## Genotype GVCFs across all samples simultaneously --
## Only analyzes variants on contigs of minimum length specified in contig filename
ref=/scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fna

while read -a line
do

if [ ! -f "final_filtered_${line[0]}_nobaseQrecal_SNPs_nearindelfilt.vcf" ]; then

gatk --java-options "-Xmx60G" GenomicsDBImport \
	-V 3811_S67_nobaseQrecal.g.vcf \
	-V 3850_S45_nobaseQrecal.g.vcf \
	-V 3910_S63_nobaseQrecal.g.vcf \
	-V 4058_S46_nobaseQrecal.g.vcf \
	-V 4100_S68_nobaseQrecal.g.vcf \
	-V 4195_S23_nobaseQrecal.g.vcf \
	-V 4444_S24_nobaseQrecal.g.vcf \
	-V 4459_S64_nobaseQrecal.g.vcf \
	-V 4793_S47_nobaseQrecal.g.vcf \
	-V 4795_S25_nobaseQrecal.g.vcf \
	-V 4796_S26_nobaseQrecal.g.vcf \
	-V 4816_S27_nobaseQrecal.g.vcf \
	-V 4835_S56_nobaseQrecal.g.vcf \
	-V 4867_S65_nobaseQrecal.g.vcf \
	-V 4890_S48_nobaseQrecal.g.vcf \
	-V 4893_S49_nobaseQrecal.g.vcf \
	-V 4895_S57_nobaseQrecal.g.vcf \
	-V 4897_S28_nobaseQrecal.g.vcf \
	-V 4901_S29_nobaseQrecal.g.vcf \
	-V 4910_S50_nobaseQrecal.g.vcf \
	-V 4915_S51_nobaseQrecal.g.vcf \
	-V 4919_S30_nobaseQrecal.g.vcf \
	-V 4922_S58_nobaseQrecal.g.vcf \
	-V 4927_S59_nobaseQrecal.g.vcf \
	-V 4936_S31_nobaseQrecal.g.vcf \
	-V 4943_S32_nobaseQrecal.g.vcf \
	-V 4946_S33_nobaseQrecal.g.vcf \
	-V 4950_S66_nobaseQrecal.g.vcf \
	-V 4952_S34_nobaseQrecal.g.vcf \
	-V 4960_S60_nobaseQrecal.g.vcf \
	-V 4962_S52_nobaseQrecal.g.vcf \
	-V 4970_S61_nobaseQrecal.g.vcf \
	-V 4976_S53_nobaseQrecal.g.vcf \
	-V 5018_S35_nobaseQrecal.g.vcf \
	-V 5026_S62_nobaseQrecal.g.vcf \
	-V 5038_S36_nobaseQrecal.g.vcf \
	-V 5039_S37_nobaseQrecal.g.vcf \
	-V 5040_S54_nobaseQrecal.g.vcf \
	-V 5046_S38_nobaseQrecal.g.vcf \
	-V 5050_S39_nobaseQrecal.g.vcf \
	-V 5054_S55_nobaseQrecal.g.vcf \
	-V 5060_S40_nobaseQrecal.g.vcf \
	-V 5075_S41_nobaseQrecal.g.vcf \
	-V 5114_S42_nobaseQrecal.g.vcf \
	-V 5123_S43_nobaseQrecal.g.vcf \
	-V 545_S21_nobaseQrecal.g.vcf \
	-V 562_S44_nobaseQrecal.g.vcf \
	-V 572_S22_nobaseQrecal.g.vcf \
	--genomicsdb-workspace-path ${line[0]}_database_gatk_genomics \
	--intervals ${line[0]}

gatk --java-options "-Xmx60G" GenotypeGVCFs \
	-R $ref \
	-V gendb://${line[0]}_database_gatk_genomics \
	-O ${line[0]}_genotype_output_nobaseQrecal.vcf
	
gatk --java-options "-Xmx60G" VariantFiltration \
	-R $ref \
	-V ${line[0]}_genotype_output_nobaseQrecal.vcf \
	-O ${line[0]}_genotype_output_nobaseQrecal_filtered.vcf \
	--filter-name "QD" \
	--filter-expression "QD < 2.0" \
	--filter-name "FS" \
	--filter-expression "FS > 60.0" \
	--filter-name "SOR" \
	--filter-expression "SOR > 5.0" \
	--filter-name "MQ" \
	--filter-expression "MQ < 40.0" \
	--filter-name "MQRankSum" \
	--filter-expression " MQRankSum < -3.0 || MQRankSum > 3.0" \
	--filter-name "ReadPosRankSum" \
	--filter-expression "ReadPosRankSum < -8.0"

echo ${line[0]} >> /_group_/_group__done_contigs.txt


## --------------------------------
## Finally, remove SNPs within 5 bp of indels using bcftools before retaining only 
## biallelic SNPs with SelectVariants

bcftools filter --SnpGap 5 \
--threads 4 \
-o ${line[0]}_genotype_output_nobaseQrecal_filtered_nearindelfilt.vcf \
${line[0]}_genotype_output_nobaseQrecal_filtered.vcf

gatk --java-options "-Xmx60G" SelectVariants \
	-R $ref \
	-V ${line[0]}_genotype_output_nobaseQrecal_filtered_nearindelfilt.vcf \
	--select-type-to-include SNP \
	-select 'vc.isNotFiltered()' \
	-restrict-alleles-to BIALLELIC \
	-O final_filtered_${line[0]}_nobaseQrecal_SNPs_nearindelfilt.vcf

echo ${line[0]} >> /_group_/_group__round2_done_contigs.txt

fi

done < /home/amh0254/krat_roh_analyses/sample_lists/_group_.txt

