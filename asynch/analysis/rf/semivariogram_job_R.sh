#!/bin/bash
# Job name:
#SBATCH --job-name=semivar
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=12:00:00
#
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Command(s) to run (example):

stdbuf -i0 -o0 -e0 command

module load r r-packages r-spatial

Rscript --vanilla ./assess_spat_autocorr.r > semivariogram.Rout