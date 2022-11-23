#!/bin/bash

#SBATCH --job-name=02__group__krat
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's partition: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_02_read_mapping

## Set batch group
GROUP=_group_

## Set reference genome
REF=/scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fna.gz

## Create a directory on /scratch
mkdir /scratch/${USER}/${PROJ}/_group_

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/_group_

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/_group_/


## --------------------------------
## Load modules
module load bwa/0.7.17
module load samtools/1.11
module load picard/2.23.9
#module load bedtools/2.29.2


## --------------------------------
## Align reads to default (smaller) genome - keep lanes separate for now
while read -a line
do

# bwa mem -t 20 -M ${REF} \
# /scratch/avrilh/kratroh_01_read_qc/_group_/${line[0]}_L002_R1_001_val_1.fq.gz \
# /scratch/avrilh/kratroh_01_read_qc/_group_/${line[0]}_L002_R2_001_val_2.fq.gz \
# > ${line[0]}_L002.bam
# 
# ## add read group information
# java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
# I=${line[0]}_L002.bam \
# O=sorted_${line[0]}_L002_rgroups.bam \
# SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=${line[0]} RGPU=unit2


## get some mapping stats
samtools flagstat -@ 20 -O tsv sorted_${line[0]}_L002_rgroups.bam > ${line[0]}_rgroups_flagstat_out.txt

done < /home/amh0254/krat_roh_analyses/sample_lists/_group_.txt


## --------------------------------

## Copy results back to project output directory (in home)
cp sorted_*_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files/
cp *_rgroups_flagstat_out.txt /home/amh0254/krat_roh_analyses/02_read_mapping/flagstat_out/

