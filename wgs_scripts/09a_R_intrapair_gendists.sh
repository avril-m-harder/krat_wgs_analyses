#!/bin/bash

#SBATCH --job-name=R_gendists
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com


##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_09_R_analyses

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
# chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/


# -----------------------------------------------------------------------------
# Load modules and copy files
# -----------------------------------------------------------------------------
module load R/4.2.1

# cp /home/amh0254/krat_roh_analyses/09_R_analyses/input/* .
cp /home/amh0254/krat_roh_analyses/scripts/R_scripts/* .


# -----------------------------------------------------------------------------
# Run R script
# -----------------------------------------------------------------------------
Rscript EASLEY_01_intrapair_gendist_calcs.R


# -----------------------------------------------------------------------------
# Copy output to /home/
# -----------------------------------------------------------------------------
# cp intrapair_gen_distances.txt /home/amh0254/krat_roh_analyses/09_R_analyses/output/

cp intrapair_gen_distances_w_sharedwindIDs.txt \
/home/amh0254/krat_roh_analyses/09_R_analyses/output/