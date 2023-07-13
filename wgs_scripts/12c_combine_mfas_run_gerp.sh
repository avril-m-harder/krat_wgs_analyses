#!/bin/bash

#SBATCH --job-name=comb_mfa
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mem 32000M
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

mkdir /scratch/${USER}/${PROJ}/mfa_prep/dspec/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/mfa_prep/


## --------------------------------
## Load modules
module load seqtk/1.3
module load singularity/3.8.4


## --------------------------------
## Loop over list of species files, add contig seq to appropriate file, and rename to 
## include species names.
## Final set of files should comprise 1 file per contig, 1 entry per species with data for
## that contig.
while read CON
do
	cd /scratch/${USER}/${PROJ}/mfa_prep/	

	## get D. spectabilis sequence and format name
	echo ${CON} > ${CON}.txt	
	echo -e "\n" >> ${CON}.txt
	
# 	seqtk subseq \
# 	/scratch/avrilh/kratroh_01_assembindex/simple_contig_names_dspec_genbank_assem.fna \
# 	${CON}.txt > ./dspec/${CON}.mfa
# 	
# 	sed -i "s/${CON}/Dipodomys_spectabilis/g" ./dspec/${CON}.mfa
	
	## get all the other species' sequences
	find . -name "${CON}.mfa" -exec cat {} + > combined_${CON}.mfa
	echo "${CON} "
	grep ">" combined_${CON}.mfa | wc -l
	
	## oops, rename sequences in .mfa files to match tree names
	sed -i "s/_${CON}//g" combined_${CON}.mfa
	
	cd ../gerp_output/
	
# 	singularity exec ../gerp_program/gerp_2.1--hfc679d8_0.sif gerpcol \
# 	-a \
# 	-v \
# 	-z \
# 	-t /home/amh0254/krat_roh_analyses/12_gerp/rodent_species_for_timetree.nwk \
# 	-f ../mfa_prep/combined_${CON}.mfa \
# 	-e Dipodomys_ordii
# 	
# 	mv ../mfa_prep/combined_${CON}.mfa.rates ./ordii_ref_combined_${CON}.mfa.rates
		
	singularity exec ../gerp_program/gerp_2.1--hfc679d8_0.sif gerpcol \
	-a \
	-v \
	-z \
	-t /home/amh0254/krat_roh_analyses/12_gerp/rodent_species_for_timetree.nwk \
	-f ../mfa_prep/combined_${CON}.mfa \
	-e Dipodomys_spectabilis
	
	mv ../mfa_prep/combined_${CON}.mfa.rates ./spectabilis_ref_combined_${CON}.mfa.rates
	
	rm ../mfa_prep/combined_${CON}.mfa
	
done < /home/amh0254/krat_roh_analyses/sample_lists/filtered_1Mb_contigs.txt
