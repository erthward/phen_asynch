#!/bin/bash
# Job name:
#SBATCH --job-name=ch3_phen_val
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

module load gdal python

python /global/home/users/drewhart/seasonality/seasonal_asynchrony/phen/validation/compare_NIRv_SIF_maps/compare_NIRv_SIF_fitted_phenology.py

