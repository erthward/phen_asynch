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
# run for both vars and for all neigh rads
for var in NIRv SIF
do
   for neigh_rad in 50 100 150
   do
      echo "Running for var ${var}, neigh_rad ${neigh_rad}"
      Rscript --vanilla /global/home/users/drewhart/seasonality/seasonal_asynchrony/analysis/asynch_corr/asynch_corr.r $var $neigh_rad > ch3_rf_${var}_${neigh_rad}.Rout 
 done
done

