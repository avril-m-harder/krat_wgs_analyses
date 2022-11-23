#!/bin/bash

#SBATCH --job-name=krat_assembindex
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_01_assembindex

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/


## --------------------------------
## Load modules
module load samtools/1.11
module load picard/2.23.9
module load bwa/0.7.17


## --------------------------------
## Download assembly as published
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/054/845/GCA_019054845.1_ASM1905484v1/GCA_019054845.1_ASM1905484v1_genomic.fna.gz
# 
# mv GCA_019054845.1_ASM1905484v1_genomic.fna.gz dspec_genbank_assem.fna.gz


## --------------------------------
## Unzip FASTA and create indexes needed for GATK Best Practices 
# bwa index dspec_genbank_assem.fna.gz
# 
# samtools faidx dspec_genbank_assem.fna.gz
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar CreateSequenceDictionary \
# REFERENCE=dspec_genbank_assem.fna \
# OUTPUT=dspec_genbank_assem.fna.dict


## --------------------------------
## Unzip FASTA and create indexes needed for GATK Best Practices
gunzip dspec_genbank_assem.fna.gz

bwa index dspec_genbank_assem.fna

samtools faidx dspec_genbank_assem.fna
 
java -jar /tools/picard-2.23.9/libs/picard.jar CreateSequenceDictionary \
REFERENCE=dspec_genbank_assem.fna \
OUTPUT=dspec_genbank_assem.fna.dict