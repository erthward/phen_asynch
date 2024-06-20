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

module load gdal/2.2.3
module load r/4.0.3
module load r-packages
module load r-spatial/2020-11-30-r40

# run for both vars and for all neigh rads
for var in NIRv SIF
do
   for neigh_rad in 50 100 150
   do
      echo
      echo "=========================================================="
      echo
      echo "Prepping data for var ${var}, neigh_rad ${neigh_rad}..."
      echo
      Rscript --vanilla /global/home/users/drewhart/seasonality/seasonal_asynchrony/asynch/anal/drivers/prep_data/prep_phen_asynch_rf_data.r $var $neigh_rad > ch3_rf_data_prep_${var}_${neigh_rad}.Rout 
      for coords_as_covars in y n
      do
        echo
        echo ".........................................................."
        echo
        echo "Running random forest model for var ${var}, neigh_rad ${neigh_rad}, coord_as_covars ${coords_as_covars}..."
        echo
         Rscript --vanilla /global/home/users/drewhart/seasonality/seasonal_asynchrony/asynch/anal/drivers/run_rf/run_phen_asynch_rf.r $var $neigh_rad $coords_as_covars > ch3_rf_${var}_${neigh_rad}_${coords_as_covars}.Rout 
      done
 done
done

