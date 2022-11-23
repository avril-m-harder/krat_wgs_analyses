#!/bin/bash

#SBATCH --job-name=plink
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 300:00:00
#SBATCH --mem=32000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set project name
PROJ=kratroh_05_plink

## Set round number
RD=3

## Create a directory on /scratch
mkdir /scratch/${USER}/${PROJ}/round${RD}/

# ## Set permissions for directory
# chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/round${RD}/

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load plink/1.9


# ----------------------------------------------------------------------
# Convert VCF to PLINK format(s)
# ----------------------------------------------------------------------

# plink \
# 	--vcf ../../kratroh_03_gatk/krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz \
# 	--allow-extra-chr \
# 	--out krat_plink


# -----------------------------------------------------------------------------
# PLINK parameters - define arrays for each of the PLINK parameters we want to
# vary and test.
# -----------------------------------------------------------------------------

declare -a phwh=(3)                       # Values for -homozyg-window-het
declare -a phwm=(2 5 10 15 20)                  # Values for -homozyg-window-missing
declare -a phws=(50)                      # Values for -homozyg-window-snp
declare -a phzd=(50)                      # Values for -homozyg-density
declare -a phzg=(1000)                    # Values for -homozyg-gap
declare -a phwt=(0.05)                    # Values for -homozyg-window-threshold
declare -a phzs=(50)                      # Values for -homozyg-snp
declare -a phzk=(100)                     # Values for -homozyg-kb

# ----------------------------------------------------------------------
# Run PLINK to ID ROHs
# ----------------------------------------------------------------------

for p1 in ${phwh[@]}; do                             # -homozyg-window-het
	for p2 in ${phwm[@]}; do                         # -homozyg-winodow-missing
		for p3 in ${phws[@]}; do                     # -homozyg-window-snp
			for p4 in ${phzd[@]}; do                 # -homozyg-density
				for p5 in ${phzg[@]}; do             # -homozyg-gap
					for p6 in ${phwt[@]}; do         # -homozyg-window threshold
						for p7 in ${phzs[@]}; do     # -homozyg-snp
							for p8 in ${phzk[@]}; do # -homozyg-kb
								PARAM_SUFFIX=_phwh_${p1}_phwm_${p2}_phws_${p3}_phzd_${p4}_phzg_${p5}_phwt_${p6}_phzs_${p7}_phzk_${p8}
								IN_FILE=krat_plink
								OUT_FILE=krat_plink_roh${PARAM_SUFFIX}

								plink --bim ${IN_FILE}.bim \
									--bed ${IN_FILE}.bed \
									--fam ${IN_FILE}.fam \
									--allow-extra-chr \
									--homozyg-window-het ${p1} \
									--homozyg-window-missing ${p2} \
									--homozyg-window-snp ${p3} \
									--homozyg-density ${p4} \
									--homozyg-gap ${p5} \
									--homozyg-window-threshold ${p6} \
									--homozyg-snp ${p7} \
									--homozyg-kb ${p8} \
									--out ${OUT_FILE}
							done
						done
					done
				done
			done
		done
	done
done

mkdir round${RD}_hom_files
mv *.hom round${RD}_hom_files