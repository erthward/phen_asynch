#!/bin/bash
# Job name:
#SBATCH --job-name=asynch
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

# loop over neigh_rads and vars, calculating asynch maps for each
#for neigh_rad in 50 100 150
for neigh_rad in 50
do
  #for var in NIRv NIRv_STRICT SIF SIF_STRICT tmmn pr def cloud 
   for var in NIRv cloud
   do
      echo "() () () () () () () () () () () () () () () () () () () () () () () () ()\n"
      echo "NOW PROCESSING ${var} FOR NEIGH RAD ${neigh_rad} METERS...\n"
      julia -p 32 /global/home/users/drewhart/seasonality/seasonal_asynchrony/asynch/calculation/calc_asynch_REDO.jl --neigh_rad $neigh_rad --var var > ch3_asynch_job.jlout
      echo "FINISHED PROCESSING ${var} FOR NEIGH RAD ${neigh_rad} METERS.\n"
   done
done
