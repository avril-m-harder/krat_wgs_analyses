#!/bin/bash

#SBATCH --job-name=snpeff_annot
#SBATCH --partition=general 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mem 32000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

### J's partition: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_10_snpeff

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/



# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load snpeff/5.0
module load bcftools/1.15
module load vcftools/0.1.17
module load perl 


# -----------------------------------------------------------------------------
# Copy full and LD-pruned VCF files, rename chromosomes to match RefSeq names
# -----------------------------------------------------------------------------
# cp /home/amh0254/krat_roh_analyses/03_gatk/krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz ./vcf_files/
# cp /home/amh0254/krat_roh_analyses/03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz ./vcf_files/
# cp /home/amh0254/krat_roh_analyses/sample_lists/genbank_refseq_chrom_names.txt \
# ./vcf_files/
# 
# bcftools annotate \
# --rename-chrs ./vcf_files/genbank_refseq_chrom_names.txt \
# -o ./vcf_files/krat_final_LDpruned_refseqchromnames.vcf.gz \
# -O z \
# ./vcf_files/krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz
# 
# bcftools annotate \
# --rename-chrs ./vcf_files/genbank_refseq_chrom_names.txt \
# -o ./vcf_files/krat_final_refseqchromnames.vcf.gz \
# -O z \
# ./vcf_files/krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz


# -----------------------------------------------------------------------------
# Run SnpEff
# -----------------------------------------------------------------------------
# java -Xmx32g -jar /tools/snpeff-5.0/snpEff.jar \
# dspec \
# ./vcf_files/krat_final_LDpruned_refseqchromnames.vcf.gz > \
# dspec_n48_snpeff_ann_LDpruned.vcf.gz
# 
# java -Xmx32g -jar /tools/snpeff-5.0/snpEff.jar \
# dspec \
# ./vcf_files/krat_final_refseqchromnames.vcf.gz > \
# dspec_n48_snpeff_ann.vcf.gz


# -----------------------------------------------------------------------------
# Filter results to allow zero missingness
# -----------------------------------------------------------------------------
# ## zero missingness
# vcftools --gzvcf dspec_n48_snpeff_ann_LDpruned.vcf.gz \
# --recode \
# --recode-INFO-all \
# --max-missing 1 \
# --out dspec_n48_snpeff_ann_LDpruned_nomissingness
# 
# vcftools --gzvcf dspec_n48_snpeff_ann.vcf.gz \
# --recode \
# --recode-INFO-all \
# --max-missing 1 \
# --out dspec_n48_snpeff_ann_nomissingness
# 
# ## split multiple annotations into single entries
# cat dspec_n48_snpeff_ann_LDpruned_nomissingness.recode.vcf | \
# /home/amh0254/krat_roh_analyses/scripts/vcfEffOnePerLine.pl \
# > dspec_n48_snpeff_ann_LDpruned_nomissingness_splitannots.recode.vcf
# 
# cat dspec_n48_snpeff_ann_nomissingness.recode.vcf | \
# perl /home/amh0254/krat_roh_analyses/scripts/vcfEffOnePerLine.pl > \
# dspec_n48_snpeff_ann_nomissingness_splitannots.recode.vcf

## use SnpSift to organize annotations
java -Xmx32g -jar /tools/snpeff-5.0/SnpSift.jar \
extractFields -s "|" -e "." \
dspec_n48_snpeff_ann_LDpruned_nomissingness_splitannots.recode.vcf \
CHROM POS REF ALT \
"ANN[*].EFFECT" "ANN[*].IMPACT:" "ANN[*].GENE:" "ANN[*].GENEID:" \
"ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE:" "ANN[*].ERRORS" > \
dspec_n48_snpeff_ann_LDpruned_nomissingness_snpsift.txt

java -Xmx32g -jar /tools/snpeff-5.0/SnpSift.jar \
extractFields -s "|" -e "." \
dspec_n48_snpeff_ann_nomissingness_splitannots.recode.vcf \
CHROM POS REF ALT \
"ANN[*].EFFECT" "ANN[*].IMPACT:" "ANN[*].GENE:" "ANN[*].GENEID:" \
"ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE:" "ANN[*].ERRORS" > \
dspec_n48_snpeff_ann_nomissingness_snpsift.txt


## compress VCFs
# gzip dspec_n48_snpeff_ann_LDpruned_nomissingness.recode.vcf
# gzip dspec_n48_snpeff_ann_LDpruned_nomissingness_splitannots.recode.vcf
# 
# gzip dspec_n48_snpeff_ann_nomissingness.recode.vcf
# gzip dspec_n48_snpeff_ann_nomissingness_splitannots.recode.vcf