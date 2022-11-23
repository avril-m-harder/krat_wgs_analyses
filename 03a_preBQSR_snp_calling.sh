#!/bin/bash

#SBATCH --job-name=_group_haplocall
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=80000
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
chmod 700 /scratch/${USER}/${PROJ}/

## cd into directory
cd /scratch/${USER}/${PROJ}/

ref=/scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fna

## --------------------------------
## Load modules
# module load samtools/1.11
module load picard/2.23.9
module load gatk/4.1.9.0


## --------------------------------
## Mark duplicate reads - GATK tools ignore them, no need to remove, just flag
while read -a line
do

# java -jar /tools/picard-2.23.9/libs/picard.jar MarkDuplicates \
# I=/scratch/avrilh/kratroh_02_read_mapping/_group_/sorted_${line[0]}_rgroups.bam \
# O=${line[0]}_small_genome_rgroups_dupmarked.bam \
# M=${line[0]}_rgroups_markdup_metrics.txt \
# MAX_RECORDS_IN_RAM=250000
# 
# ## Fix mate information
# java -jar /tools/picard-2.23.9/libs/picard.jar FixMateInformation \
# I=${line[0]}_small_genome_rgroups_dupmarked.bam \
# O=${line[0]}_small_genome_rgroups_dupmarked_fixmate.bam
# 
## Index BAM files
if [ ! -f "${line[0]}_nobaseQrecal.g.vcf.idx" ]; then

	java -jar /tools/picard-2.23.9/libs/picard.jar BuildBamIndex \
	I=${line[0]}_small_genome_rgroups_dupmarked_fixmate.bam

## Run HaplotypeCaller in GVCF mode
	gatk HaplotypeCaller \
		-R $ref \
		-I ${line[0]}_small_genome_rgroups_dupmarked_fixmate.bam \
		-stand-call-conf 20.0 \
		--emit-ref-confidence GVCF \
		-O ${line[0]}_nobaseQrecal.g.vcf

fi
	
done < /home/amh0254/krat_roh_analyses/sample_lists/_group_.txt
