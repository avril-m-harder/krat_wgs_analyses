#!/bin/bash

#SBATCH --job-name=BWA_TEST_methylseq
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=80000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set dir name
DIR=ms_kratroh_02_nfcore_methylseq_TEST

## Set project names
PROJ1=bismark
PROJ2=bwameth


## --------------------------------
## Activate conda environment and load modules

module load python/anaconda/3.10.9

source activate nextflow
module load java/15.0.1
module load singularity/3.8.4


## --------------------------------
## Run nf-score methylseq pipeline with bwa-meth

mkdir /scratch/${USER}/${DIR}/${PROJ2}/
cd /scratch/${USER}/${DIR}/${PROJ2}/

mkdir results
mkdir data
cp /scratch/${USER}/${DIR}/${PROJ1}/data/* ./data/

nextflow run nf-core/methylseq \
--input data/krat_ms_samplesheet_TEST.csv \
--outdir results \
--max_cpus 2 \
--aligner bwameth \
--fasta /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa \
--save_reference \
-resume \
-profile singularity

conda deactivate
