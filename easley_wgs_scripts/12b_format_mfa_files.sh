#!/bin/bash

#SBATCH --job-name=mfa_NUM
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
PROJ=kratroh_12_gerp

## Create a directory on /scratch
mkdir /scratch/${USER}/${PROJ}/

mkdir /scratch/${USER}/${PROJ}/mfa_prep

mkdir /scratch/${USER}/${PROJ}/mfa_prep/NUM

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/mfa_prep/NUM

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/mfa_prep/NUM
mkdir mfa_files


## --------------------------------
## Load modules
module load seqtk/1.3


## --------------------------------
## Loop over list of species files, add contig seq to appropriate file, and rename to 
## include species names.
## Final set of files should comprise 1 file per contig, 1 entry per species with data for
## that contig.
while read SP CON IT
do
	echo ${CON} > list.txt
	
	echo -e "\n" >> list.txt
	seqtk subseq \
	/scratch/avrilh/kratroh_12_gerp/consensus_fasta_files/${SP}_consensus.fna \
	list.txt > temp.fa
	
	sed -i "s/${CON}/${SP}_${CON}/g" temp.fa
	
	cat temp.fa >> ${CON}.mfa
	
	echo "${IT}"
	
done < /home/amh0254/krat_roh_analyses/12_gerp/acc_contig_lists/acc_contig_list_NUM.txt
