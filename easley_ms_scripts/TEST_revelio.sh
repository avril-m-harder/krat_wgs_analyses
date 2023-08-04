#!/bin/bash

#SBATCH --job-name=TEST_revelio
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 300:00:00
#SBATCH --mem=16000
#SBATCH --output=TEST_revelio-%j.out
#SBATCH --error=error-TEST_revelio-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set project name
PROJ=TEST_ms_kratroh_03_snp_calling

mkdir /scratch/${USER}/${PROJ}/
chmod 700 /scratch/${USER}/${PROJ}


# -----------------------------------------------------------------------------
# Activate conda environment and load modules
# -----------------------------------------------------------------------------

module load python/anaconda/3.10.9
# export CONDA_PKGS_DIRS=~/.conda/pkgs
# conda create -n revelio
source activate revelio
conda install -c bioconda pysam

source /home/amh0254/krat_roh_analyses/scripts/easley_ms_scripts/logging_functions.sh


# -----------------------------------------------------------------------------
# Copy test data and subset to WGS-filtered contigs
# -----------------------------------------------------------------------------

cd /scratch/${USER}/${PROJ}

module load samtools/1.17

## *.sorted.bam files have been through markdups and sort
cp /scratch/avrilh/ms_kratroh_02_nfcore_methylseq/bismark/work/86/\
d3b8650c59b4a96bae568d73cc50fc/1568.sorted.bam .

samtools view \
--threads 3 \
-L /home/amh0254/krat_roh_analyses/sample_lists/filtered_contigs.bed \
-b \
1568.sorted.bam \
> 1568.sorted.filtered.bam


# -----------------------------------------------------------------------------
# Create masked BAM file with Revelio
# -----------------------------------------------------------------------------

start_logging "Run Revelio - /scratch/${USER}/${PROJ}/"

samtools calmd \
--threads 3 \
-b 1568.sorted.filtered.bam \
/scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa 1> \
1568.sorted.filtered.calmd.bam \
2> /dev/null

/home/amh0254/revelio-main/revelio.py \
1568.sorted.filtered.calmd.bam \
1568.sorted.filtered.masked.bam

# -----------------------------------------------------------------------------
# [Add in GATK steps that match(?) what I did for WGS samples]
# -----------------------------------------------------------------------------



stop_logging

conda deactivate