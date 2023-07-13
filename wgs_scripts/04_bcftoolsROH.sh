#!/bin/bash

#SBATCH --job-name=krat_roh
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=40000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_04_bcftoolsroh

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/
# 
# ## Set permissions for directory
# chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/


## --------------------------------
## Load modules
module load bcftools/1.11
# module load samtools/1.11
module load vcftools/0.1.17
module load htslib/1.11


## --------------------------------
## Generate allele frequency file for use with bcftools/ROH & index

bgzip -@ 20 -c \
../kratroh_03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf \
> ../kratroh_03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
../kratroh_03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz \
| bgzip -c > freqs_krat_final_allfiltcontigs_all_samps.tab.gz

tabix -s1 -b2 -e2 freqs_krat_final_allfiltcontigs_all_samps.tab.gz


## --------------------------------
## Identify ROHs using both GT and PL settings

bcftools roh \
--threads 20 \
-o krat_final_allfiltcontigs_all_samps_PL.txt \
--AF-file freqs_krat_final_allfiltcontigs_all_samps.tab.gz \
../kratroh_03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz

bcftools roh \
--threads 20 \
--GTs-only 30 \
-o krat_final_allfiltcontigs_all_samps_GT.txt \
--AF-file freqs_krat_final_allfiltcontigs_all_samps.tab.gz \
../kratroh_03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz

## --------------------------------
## Extract information on ROHs (i.e., exclude information on individual sites)

grep "^RG" krat_final_allfiltcontigs_all_samps_GT.txt > \
krat_final_allfiltcontigs_all_samps_GT_RG_ONLY.txt

grep "^RG" krat_final_allfiltcontigs_all_samps_PL.txt > \
krat_final_allfiltcontigs_all_samps_PL_RG_ONLY.txt





## --------------------------------
## bcftools roh details
## other info: https://samtools.github.io/bcftools/howtos/roh-calling.html
##
## About:   HMM model for detecting runs of autozygosity.
## Usage:   bcftools roh [options] <in.vcf.gz>
## 
## General Options:
##         --AF-dflt <float>              if AF is not known, use this allele frequency [skip]
##         --AF-tag <TAG>                 use TAG for allele frequency
##         --AF-file <file>               read allele frequencies from file (CHR\tPOS\tREF\tALT\tAF)
##     -b  --buffer-size <int[,int]>      buffer size and the number of overlapping sites, 0 for unlimited [0]
##                                            If the first number is negative, it is interpreted as the maximum memory to
##                                            use, in MB. The default overlap is set to roughly 1% of the buffer size.
##     -e, --estimate-AF [TAG],<file>     estimate AF from FORMAT/TAG (GT or PL) of all samples ("-") or samples listed
##                                             in <file>. If TAG is not given, the frequency is estimated from GT by default
##         --exclude <expr>               exclude sites for which the expression is true
##     -G, --GTs-only <float>             use GTs and ignore PLs, instead using <float> for PL of the two least likely genotypes.
##                                            Safe value to use is 30 to account for GT errors.
##         --include <expr>               select sites for which the expression is true
##     -i, --ignore-homref                skip hom-ref genotypes (0/0)
##     -I, --skip-indels                  skip indels as their genotypes are enriched for errors
##     -m, --genetic-map <file>           genetic map in IMPUTE2 format, single file or mask, where string "{CHROM}"
##                                            is replaced with chromosome name
##     -M, --rec-rate <float>             constant recombination rate per bp
##     -o, --output <file>                write output to a file [standard output]
##     -O, --output-type [srz]            output s:per-site, r:regions, z:compressed [sr]
##     -r, --regions <region>             restrict to comma-separated list of regions
##     -R, --regions-file <file>          restrict to regions listed in a file
##     -s, --samples <list>               list of samples to analyze [all samples]
##     -S, --samples-file <file>          file of samples to analyze [all samples]
##     -t, --targets <region>             similar to -r but streams rather than index-jumps
##     -T, --targets-file <file>          similar to -R but streams rather than index-jumps
##         --threads <int>                use multithreading with <int> worker threads [0]
## 
## HMM Options:
##     -a, --hw-to-az <float>             P(AZ|HW) transition probability from HW (Hardy-Weinberg) to AZ (autozygous) state [6.7e-8]
##     -H, --az-to-hw <float>             P(HW|AZ) transition probability from AZ to HW state [5e-9]
##     -V, --viterbi-training <float>     estimate HMM parameters, <float> is the convergence threshold, e.g. 1e-10 (experiment