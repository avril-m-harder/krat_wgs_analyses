#!/bin/bash

#SBATCH --job-name=gerp_dl
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_11_gerp

## Create a directory on /scratch
mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/

## --------------------------------
## Load module, create, and activate conda environment
module load python/anaconda/3.8.6
# 
# export CONDA_PKGS_DIRS=~/.conda/pkgs
# 
# conda create -n gerp
# 
source activate gerp
# 
# conda install -c conda-forge ncbi-datasets-cli
# 
# conda deactivate


## --------------------------------
## Download reference genomes using ncbi datasets