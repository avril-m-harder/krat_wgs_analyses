#!/bin/bash

#SBATCH --job-name=02c_krat
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's partition: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_02_read_mapping


## --------------------------------
## Load modules
# module load bwa/0.7.17
# module load samtools/1.11
module load picard/2.23.9
# module load bedtools/2.29.2


## --------------------------------
## Merge BAM files for Duke and UofI reseq
# cd /scratch/avrilh/kratroh_02_read_mapping/

## merge BAM files within each sample
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./18/sorted_5050_TTAACCTTCG-TAATGGCAAG_L002_rgroups.bam \
# I=./14/sorted_5050_S39_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_5050_S39_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./18/sorted_5060_GCGTGCTGTG-CGGTGACACC_L002_rgroups.bam \
# I=./14/sorted_5060_S40_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_5060_S40_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./18/sorted_5075_CGAGAGGCGT-GAGACATAAT_L002_rgroups.bam \
# I=./15/sorted_5075_S41_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_5075_S41_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./17/sorted_4960_CAAGTTATTG-GCGCAGAGTA_L002_rgroups.bam \
# I=./10/sorted_4960_S60_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_4960_S60_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./17/sorted_4897_AATGCGAACA-ATTACTCACC_L002_rgroups.bam \
# I=./06/sorted_4897_S28_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_4897_S28_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./17/sorted_5018_CGCTGTCTCA-ATGTCGTGGT_L002_rgroups.bam \
# I=./12/sorted_5018_S35_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_5018_S35_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./17/sorted_5038_GCTAATAGGA-AGAGCACTAG_L002_rgroups.bam \
# I=./12/sorted_5038_S36_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_5038_S36_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./18/sorted_545_TGGTAGAGAT-TGTTGTTCGT_L002_rgroups.bam \
# I=./16/sorted_545_S21_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_545_S21_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./18/sorted_572_AGTACTCATG-GTAGAGTCAG_L002_rgroups.bam \
# I=./16/sorted_572_S22_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_572_S22_rgroups.bam
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar MergeSamFiles \
# I=./17/sorted_4927_ACCAAGCAGG-TGGCTCGCAG_L002_rgroups.bam \
# I=./08/sorted_4927_S59_rgroups.bam \
# SORT_ORDER=coordinate \
# O=./combined_bams/sorted_4927_S59_rgroups.bam


## --------------------------------
## Confirm read #'s maintained in merge
# cd ./combined_bams/
# 
# while read -a line
# do
# 
# ## get some mapping stats
# samtools flagstat -@ 20 -O tsv sorted_${line[0]}_rgroups.bam > ${line[0]}_rgroups_flagstat_out.txt
# 
# done < /home/amh0254/krat_roh_analyses/sample_lists/double_seqed_samps.txt


## --------------------------------
## If #'s check out, overwrite with new merged BAM files
## and only use groups 01-16 going forward
# cd ./combined_bams/

# cp sorted_5050_S39_rgroups.bam ../14/
# cp sorted_5060_S40_rgroups.bam ../14/
# cp sorted_5075_S41_rgroups.bam ../15/
# cp sorted_4960_S60_rgroups.bam ../10/
# cp sorted_4897_S28_rgroups.bam ../06/
# cp sorted_5018_S35_rgroups.bam ../12/
# cp sorted_5038_S36_rgroups.bam ../12/
# cp sorted_545_S21_rgroups.bam ../16/
# cp sorted_572_S22_rgroups.bam ../16/
# cp sorted_4927_S59_rgroups.bam ../08/
# 
# cp sorted_5050_S39_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_5060_S40_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_5075_S41_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_4960_S60_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_4897_S28_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_5018_S35_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_5038_S36_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_545_S21_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_572_S22_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files
# cp sorted_4927_S59_rgroups.bam /home/amh0254/krat_roh_analyses/02_read_mapping/sorted_bam_files

## --------------------------------
## Added on May 12, 2022 because these sets of reads have different sample names in the 
## read groups, making them look like multisample BAM files to GATK
cd /scratch/avrilh/kratroh_03_gatk/

mv 5050_S39_small_genome_rgroups_dupmarked_fixmate.bam \
5050_S39_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=5050_S39_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=5050_S39_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=5050_S39 RGPU=unit2

mv 5060_S40_small_genome_rgroups_dupmarked_fixmate.bam \
5060_S40_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=5060_S40_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=5060_S40_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=5060_S40 RGPU=unit2

mv 5075_S41_small_genome_rgroups_dupmarked_fixmate.bam \
5075_S41_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=5075_S41_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=5075_S41_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=5075_S41 RGPU=unit2

mv 4960_S60_small_genome_rgroups_dupmarked_fixmate.bam \
4960_S60_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=4960_S60_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=4960_S60_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=4960_S60 RGPU=unit2

mv 4897_S28_small_genome_rgroups_dupmarked_fixmate.bam \
4897_S28_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=4897_S28_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=4897_S28_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=4897_S28 RGPU=unit2

mv 5018_S35_small_genome_rgroups_dupmarked_fixmate.bam \
5018_S35_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=5018_S35_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=5018_S35_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=5018_S35 RGPU=unit2

mv 5038_S36_small_genome_rgroups_dupmarked_fixmate.bam \
5038_S36_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=5038_S36_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=5038_S36_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=5038_S36 RGPU=unit2

mv 545_S21_small_genome_rgroups_dupmarked_fixmate.bam \
545_S21_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=545_S21_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=545_S21_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=545_S21 RGPU=unit2

mv 572_S22_small_genome_rgroups_dupmarked_fixmate.bam \
572_S22_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=572_S22_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=572_S22_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=572_S22 RGPU=unit2

mv 4927_S59_small_genome_rgroups_dupmarked_fixmate.bam \
4927_S59_small_genome_rgroups_dupmarked_fixmate_OLD.bam

java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=4927_S59_small_genome_rgroups_dupmarked_fixmate_OLD.bam \
O=4927_S59_small_genome_rgroups_dupmarked_fixmate.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=4927_S59 RGPU=unit2
