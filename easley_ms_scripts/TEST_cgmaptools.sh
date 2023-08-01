#!/bin/bash

#SBATCH --job-name=cgmt_TEST
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mem=16000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set project name
PROJ=TEST_ms_kratroh_03_snp_calling

mkdir /scratch/${USER}/${PROJ}/
chmod 700 /scratch/${USER}/${PROJ}
cd /scratch/${USER}/${PROJ}


# -----------------------------------------------------------------------------
# Copy test data and subset to WGS-filtered contigs
# -----------------------------------------------------------------------------

# module load samtools/1.17
# 
# cp /scratch/avrilh/ms_kratroh_02_nfcore_methylseq/bismark/work/86/\
# d3b8650c59b4a96bae568d73cc50fc/1568.sorted.bam .
# 
# samtools view \
# -L /home/amh0254/krat_roh_analyses/sample_lists/filtered_contigs.bed \
# -b \
# 1568.sorted.bam \
# > 1568.sorted.filtered.bam


# -----------------------------------------------------------------------------
# Create ATCGmap file and run SNP calling (Bayesian)
# -----------------------------------------------------------------------------

cgmaptools convert bam2cgmap \
-b 1568.sorted.bam \
-g /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa \
-o TEST_1568

# cgmaptools snv \
# -m bayes --bayes-dynamicP \
# --rmOverlap \
# -i <ATCGmap> \
# -v <VCF>


## There are two strategies we can use to call SNPs, binomial or bayesian mode. The 
## default choice is the binomial mode as it’s faster, but the overall prediction
## precision of bayesian mode is better. If you are not in hurry, we recommand you
## the beyasian mode for more accuracy callings.
## 
## Also to increase prediction precision, --bayes-dynamicP is recommanded when using 
## bayesian mode, which dynamicsly adjust p-value threshoulds according to the read 
## coverage.
