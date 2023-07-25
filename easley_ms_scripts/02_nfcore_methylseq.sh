#!/bin/bash

#SBATCH --job-name=methylseq
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set dir name
DIR=ms_kratroh_02_nfcore_methylseq

## Set project names
PROJ1=bismark


## --------------------------------
## Activate conda environment and load modules

module load python/anaconda/3.10.9

source activate nextflow
module load java/15.0.1
module load singularity/3.8.4


## --------------------------------
## Download test data

mkdir /scratch/${USER}/${DIR}/
chmod -R 700 mkdir /scratch/${USER}/${DIR}/
mkdir /scratch/${USER}/${DIR}/${PROJ1}/
cd /scratch/${USER}/${DIR}/${PROJ1}/

mkdir data
# cp /home/amh0254/krat_roh_analyses/sample_lists/krat_ms_samplesheet.csv ./data/


## --------------------------------
## Run nf-score methylseq pipeline

# mkdir results

## aligning with bismark -- requires a manual fix after this step
# nextflow run nf-core/methylseq \
# --input data/krat_ms_samplesheet.csv \
# --outdir results \
# --max_cpus 8 \
# --fasta /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa \
# --save_reference \
# -resume \
# -profile singularity

## for some reason, ref .fa can't be used when it's just a symbolic link? idk. this fixes
## it tho, alongside specif.
# cd /scratch/avrilh/ms_kratroh_02_nfcore_methylseq/bismark/work/78/9f96d97b18f0ecb078af557a878f2e/BismarkIndex/
# rm dspec_genbank_assem.fa
# cp /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa .


nextflow run nf-core/methylseq \
--input data/krat_ms_samplesheet.csv \
--outdir results \
--fasta /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa \
--save_reference \
-resume \
-c /home/amh0254/krat_roh_analyses/scripts/files/easleyconfig_amh.conf \
-profile singularity


conda deactivate
