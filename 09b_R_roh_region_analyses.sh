#!/bin/bash

#SBATCH --job-name=R_roh_regions
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8000
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

cp /home/amh0254/krat_roh_analyses/09_R_analyses/input/* .
cp /home/amh0254/krat_roh_analyses/scripts/R_scripts/* .


# -----------------------------------------------------------------------------
# Run R script
# -----------------------------------------------------------------------------
Rscript EASLEY_02_lo_and_hi_roh_region_ID.R


# -----------------------------------------------------------------------------
# Copy output to /home/
# -----------------------------------------------------------------------------
# cp intrapair_gen_distances.txt /home/amh0254/krat_roh_analyses/09_R_analyses/output/

cp *randomization_stats*.csv /home/amh0254/krat_roh_analyses/09_R_analyses/output/
