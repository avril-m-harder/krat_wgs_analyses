#!/bin/bash

##  Set username
USER=avrilh

## Set dir name
DIR=ms_kratroh_02_nfcore_methylseq

## Set project names
PROJ1=bismark


## --------------------------------
## Activate conda environment and load modules
## (For some reason, I get java-related errors if I rely on the config file for 
## these steps)

module load python/anaconda/3.10.9

source activate nextflow
module load java/15.0.1
module load singularity/3.8.4


## --------------------------------
## Download test data

# mkdir /scratch/${USER}/${DIR}/
# chmod -R 700 mkdir /scratch/${USER}/${DIR}/
# mkdir /scratch/${USER}/${DIR}/${PROJ1}/
cd /scratch/${USER}/${DIR}/${PROJ1}/

# mkdir data
# cp /home/amh0254/krat_roh_analyses/sample_lists/krat_ms_samplesheet.csv ./data/
# cp /home/amh0254/krat_roh_analyses/scripts/files/easleyconfig_amh.conf .


## --------------------------------
## Run nf-score methylseq pipeline

# mkdir results

## aligning with bismark -- requires a manual fix after this step
# nextflow run nf-core/methylseq \
# --input data/krat_ms_samplesheet.csv \
# --outdir results \
# --fasta /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa \
# --save_reference \
# -resume \
# -c easleyconfig_amh.conf \
# -profile singularity

## for some reason, ref .fa can't be used when it's just a symbolic link? idk. this fixes
## it tho, alongside specif.
# cd /scratch/avrilh/ms_kratroh_02_nfcore_methylseq/bismark/work/fc/e4a108a5ce95a58106631b57f1dd5b/BismarkIndex/
# rm dspec_genbank_assem.fa
# cp /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa .


nextflow run nf-core/methylseq \
--input data/krat_ms_samplesheet.csv \
--outdir results \
--fasta /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa \
--save_reference \
-resume \
-c easleyconfig_amh.conf \
-profile singularity


## below should be covered by config file
# conda deactivate
