#!/bin/bash

#SBATCH --job-name=dl
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mem=6000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set project name
PROJ=ms_kratroh_01_read_qc

# mkdir /scratch/${USER}/${PROJ}/
chmod 700 /scratch/${USER}/${PROJ}/
cd /scratch/${USER}/${PROJ}/
# mkdir raw_methylseq_data
cd raw_methylseq_data

# curl -o GMCF_1530_merged_tarball_1_20230718.tar.bz2 -L "https://posting.biotech.illinois.edu/posting/RushUniv/GMCF_1530_merged_tarball_1_20230718.tar.bz2?x-email=Kevin_Kunstman%40rush.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20230718%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230718T224048Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=4fbc8b0005a79b946f465efab314dcfc6fc6f4f82aafac2be42f856b9c6fec1f"
# 
# curl -o GMCF_1530_merged_tarball_2_20230718.tar.bz2 -L "https://posting.biotech.illinois.edu/posting/RushUniv/GMCF_1530_merged_tarball_2_20230718.tar.bz2?x-email=Kevin_Kunstman%40rush.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20230718%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230718T224048Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=d266058ab262a186d03362b500c1960bdbc78a829795fd5da3bc53c84ddccc7f"
# 
# curl -o GMCF_1530_merged_tarball_3_20230718.tar.bz2 -L "https://posting.biotech.illinois.edu/posting/RushUniv/GMCF_1530_merged_tarball_3_20230718.tar.bz2?x-email=Kevin_Kunstman%40rush.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20230718%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230718T224048Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=35bf2fbfe4ae18201763cf72efc2a879adeeaa526d1c6355e7f19cf856ed529f"
# 
# curl -o GMCF_1530_merged_tarball_4_20230718.tar.bz2 -L "https://posting.biotech.illinois.edu/posting/RushUniv/GMCF_1530_merged_tarball_4_20230718.tar.bz2?x-email=Kevin_Kunstman%40rush.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20230718%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230718T224048Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=c01c986ed1dd29cc86032fd994f978dea6a80e9507e67d67c7f96fa9aed47c5a"
# 
# cp /home/amh0254/krat_roh_analyses/raw_methylseq_data/checklist.txt .
# md5sum -c checklist.txt

tar -xvf GMCF_1530_merged_tarball_1_20230718.tar.bz2
tar -xvf GMCF_1530_merged_tarball_2_20230718.tar.bz2
tar -xvf GMCF_1530_merged_tarball_3_20230718.tar.bz2
tar -xvf GMCF_1530_merged_tarball_4_20230718.tar.bz2
