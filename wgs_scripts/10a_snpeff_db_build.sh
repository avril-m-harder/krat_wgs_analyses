#!/bin/bash

#SBATCH --job-name=snpeff_db
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
# If purged, set up SnpEff directories
# -----------------------------------------------------------------------------
# mkdir data
# cd data
# mkdir dspec
# cd ../

## copy config file from home (already edited)
cp /home/amh0254/krat_roh_analyses/10_snpeff/snpEff.config .

## copy and unzip genome FASTA file
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/054/845/GCF_019054845.1_ASM1905484v1/GCF_019054845.1_ASM1905484v1_genomic.fna.gz
# gunzip GCF_019054845.1_ASM1905484v1_genomic.fna.gz
# mv GCF_019054845.1_ASM1905484v1_genomic.fna ./data/dspec/sequences.fa
# 
# ## download and unzip GTF annotation file
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/054/845/GCF_019054845.1_ASM1905484v1/GCF_019054845.1_ASM1905484v1_genomic.gtf.gz
# gunzip GCF_019054845.1_ASM1905484v1_genomic.gtf.gz
# mv GCF_019054845.1_ASM1905484v1_genomic.gtf ./data/dspec/genes.gtf
# 

###### not including these data because database will not build
## get CDS FASTA file
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/054/845/GCF_019054845.1_ASM1905484v1/GCF_019054845.1_ASM1905484v1_cds_from_genomic.fna.gz
# gunzip GCF_019054845.1_ASM1905484v1_cds_from_genomic.fna.gz
# mv GCF_019054845.1_ASM1905484v1_cds_from_genomic.fna ./data/dspec/cds.fa 
# 
# ## not doing this for now because it doesn't seem to help anything ¯\_(ツ)_/¯
# mv GCF_019054845.1_ASM1905484v1_cds_from_genomic.fna ./data/dspec/temp_cds.fa
# awk -F '_' '/^>/ {$1=">"$4"_"$5 } 1' ./data/dspec/temp_cds.fa > \
# ./data/dspec/temp_cds_b.fa
# awk -F " " '{ print $1 }' ./data/dspec/temp_cds_b.fa > ./data/dspec/cds.fa
# 
# ## get protein FASTA file
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/054/845/GCF_019054845.1_ASM1905484v1/GCF_019054845.1_ASM1905484v1_protein.faa.gz
# gunzip GCF_019054845.1_ASM1905484v1_protein.faa.gz
# mv GCF_019054845.1_ASM1905484v1_protein.faa ./data/dspec/protein.fa


# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load snpeff/5.0


# -----------------------------------------------------------------------------
# Build SnpEff database
# -----------------------------------------------------------------------------

java -jar /tools/snpeff-5.0/snpEff.jar build \
-c snpEff.config \
-gtf22 \
-v \
dspec