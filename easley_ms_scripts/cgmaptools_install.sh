#!/bin/bash

#SBATCH --job-name=cgmt_install
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com


#### Install program in home directory
cd
wget https://github.com/guoweilong/cgmaptools/archive/refs/tags/v0.1.2.tar.gz
mv v0.1.2.tar.gz cgmaptools_v0.1.2.tar.gz
tar -zxf cgmaptools_v0.1.2.tar.gz
cd cgmaptools-0.1.2
sh install.sh

## added below line to .bashrc file:
# export PATH=/home/amh0254/cgmaptools-0.1.2:$PATH