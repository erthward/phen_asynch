#!/bin/bash
# Job name:
#SBATCH --job-name=asy100
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
## Commands to run:

stdbuf -i0 -o0 -e0 command

module load python julia/1.4.1

# loop over vars, calculating asynch maps for each
for neigh_rad in 100
do
   for var in NIRv NIRv_STRICT SIF SIF_STRICT tmmn tmmx pr def cloud 
   do
      echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      echo "NOW PROCESSING ${var} FOR NEIGH RAD ${neigh_rad} METERS..."
      julia -p 32 /global/home/users/drewhart/seasonality/seasonal_asynchrony/asynch/calculation/calc_asynch.jl --neigh_rad $neigh_rad --var $var > ch3_asynch_job_${var}_${neigh_rad}.jlout
      echo "FINISHED PROCESSING ${var} FOR NEIGH RAD ${neigh_rad} METERS."
   done
done
