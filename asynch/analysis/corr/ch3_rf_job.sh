#!/bin/bash
# Job name:
#SBATCH --job-name=ch3_rf
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Command(s) to run:

module load gdal r r-packages r-spatial
Rscript --vanilla /global/home/users/drewhart/seasonality/seasonal_asynchrony/analysis/asynch_corr/asynch_corr.r > ch3_rf.Rout 

