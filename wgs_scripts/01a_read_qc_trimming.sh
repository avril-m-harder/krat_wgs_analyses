#!/bin/bash

#SBATCH --job-name=01a__group__krat
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_01_read_qc

## Set batch group
GROUP=_group_

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/_group_

##  Move zipped fastq files to  dir to process
# while read -a line
# do
# 	cp /home/amh0254/krat_roh_analyses/raw_reseq_data/fastq_files/${line[0]}*.fastq.gz \
# 	/scratch/${USER}/${PROJ}/_group_/
# done < /home/amh0254/krat_roh_analyses/sample_lists/_group_.txt

## --------------------------------
## Activate conda environment and load other modules 
cd
source activate env1
module load python/3.8.6
module load fastqc/0.11.9

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/_group_/

## Run fastqc on all 4 read files per sample
# mkdir pretrim_output/
# mkdir posttrim_output/
# mkdir stdouts/
# 
# while read -a line
# do
# 	fastqc -t 1 -o pretrim_output/ \
# 	${line[0]}_L003_R1_001.fastq.gz\
# 	${line[0]}_L003_R2_001.fastq.gz\
# 	${line[0]}_L004_R1_001.fastq.gz\
# 	${line[0]}_L004_R2_001.fastq.gz
# done < /home/amh0254/krat_roh_analyses/sample_lists/_group_.txt
# 
# cp pretrim_output/* /home/amh0254/krat_roh_analyses/01_read_qc_trimming/pretrim_fastqc/

## Run TrimGalore to clean and trim adapters from reads,
## then run FastQC on the trimmed files
while read -a line
do

	/home/amh0254/programs/trimgalore/TrimGalore-0.6.6/trim_galore \
	--paired --quality 20 --phred33 \
	--fastqc_args "--outdir posttrim_output/" \
	--length 30 \
	${line[0]}_L003_R1_001.fastq.gz \
	${line[0]}_L003_R2_001.fastq.gz
	
	/home/amh0254/programs/trimgalore/TrimGalore-0.6.6/trim_galore \
	--paired --quality 20 --phred33 \
	--fastqc_args "--outdir posttrim_output/" \
	--length 30 \
	${line[0]}_L004_R1_001.fastq.gz \
	${line[0]}_L004_R2_001.fastq.gz

done < /home/amh0254/krat_roh_analyses/sample_lists/_group_.txt

cp posttrim_output/* /home/amh0254/krat_roh_analyses/01_read_qc_trimming/posttrim_fastqc/

## Deactivate conda environment
conda deactivate 

## Clean up
# mv *.out stdouts/
