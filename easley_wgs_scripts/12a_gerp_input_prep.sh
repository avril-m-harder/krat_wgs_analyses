#!/bin/bash

#SBATCH --job-name=gerp_NUM
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem 16000M
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
mkdir /scratch/${USER}/${PROJ}/consensus_fasta_files
mkdir /scratch/${USER}/${PROJ}/contig_covg_files

## Create batch directory on /scratch
mkdir /scratch/${USER}/${PROJ}/NUM

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/NUM

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/NUM

## make fastq and bam directories
mkdir fastq_files
mkdir bam_files
mkdir consensus_files


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
# conda install -c bioconda htsbox 


## --------------------------------
## Load other modules
module load bbmap/38.92
module load bwa/2.0
module load samtools/1.15
module load seqkit/0.14.0


## --------------------------------
## Loop over file of accession numbers and species names
while read ACC SP
	do
	cd /scratch/${USER}/${PROJ}/NUM/ 
	if [ -a ./consensus_files/${SP}_consensus.fna ];
	then
		echo "${SP} already processed"
	else
		## --------------------------------
		## Download and rename reference genome, split into fastq reads
		
	 datasets download genome accession ${ACC} --filename ${SP}.zip

		unzip ${SP}.zip
	
		mv ncbi_dataset/ ${ACC}/

		cd ${ACC}/data/${ACC}/

		mv *.fna ${ACC}.fna

		java -ea -Xms16g -cp /tools/bbmap-38.92/current/ jgi.ReformatReads \
		in=${ACC}.fna \
		out=/scratch/${USER}/${PROJ}/NUM/fastq_files/${ACC}.fq \
		qfake=35 \
		fastareadlen=50

 		## check if file exists and clean up dataset download if so
		cd /scratch/${USER}/${PROJ}/NUM/fastq_files/

		if [ -a ${ACC}.fq ];
		then
			rm -rf /scratch/${USER}/${PROJ}/NUM/${ACC}/
			rm /scratch/${USER}/${PROJ}/NUM/${SP}.zip
			rm /scratch/${USER}/${PROJ}/NUM/README.md
			echo "${SP} NCBI files cleaned up"
		else
			echo "${ACC}.fq did not generate"
		fi

		## --------------------------------
		## Align generated reads to D. spectabilis reference

		bwa mem -B 3 \
		-O 4,4 \
		-t 5 \
		/scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fna.gz \
		${ACC}.fq \
		> ../bam_files/${ACC}.bam

		cd ../bam_files/

		samtools stats ${ACC}.bam > ${ACC}_stats.txt

		## get rid of reads that did not map or alignments that are secondary
		samtools view -b -F 0x104 \
		--threads 4 \
		--output filt_${ACC}.bam \
		${ACC}.bam

		samtools sort \
		-@ 4 \
		-o sort_filt_${ACC}.bam \
		filt_${ACC}.bam

		samtools stats sort_filt_${ACC}.bam > sort_filt_${ACC}_stats.txt


		## --------------------------------
		## Convert bam to fasta
		htsbox pileup -R -q 30 -Q 30 -l 35 -s 1 \
		-b /home/amh0254/krat_roh_analyses/12_gerp/contigs_1e6.bed \
		-f /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fna.gz \
		sort_filt_${ACC}.bam > \
		../consensus_files/${SP}_consensus.fna
		
		seqkit fx2tab --length --name --header-line \
		../consensus_files/${SP}_consensus.fna > \
		../consensus_files/${SP}_contigcov.txt

		## clean up files if successful
# 		if [ -a ../consensus_files/${SP}_consensus.fna ];
# 		then
# 			rm *${ACC}.bam
# 			rm ../fastq_files/${ACC}.fq
# 			echo "${SP} BAM files cleaned up"
# 		else
# 			echo "${SP}_consensus.fna did not generate"
# 		fi
	fi
	
	cp ../consensus_files/*_contigcov.txt \
	/scratch/avrilh/kratroh_12_gerp/contig_covg_files
	
	cp ../consensus_files/*.fna \
	/scratch/avrilh/kratroh_12_gerp/consensus_fasta_files

done < /home/amh0254/krat_roh_analyses/sample_lists/rodentia_acc_spp_NUM.txt


conda deactivate