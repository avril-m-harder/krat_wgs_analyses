#!/bin/bash

#SBATCH --job-name=reformat
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

cd /scratch/avrilh/kratroh_01_assembindex

awk '/^>/ { print $1; next; }; { print; }' dspec_genbank_assem.fna > \
simple_contig_names_dspec_genbank_assem.fna
