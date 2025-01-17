#!/bin/bash
# Job name:
#SBATCH --job-name=asy50
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

# loop over vars, calculating asynch maps for each
for neigh_rad in 50
do
   for var in NIRv NIRv_STRICT SIF SIF_STRICT tmmn tmmx pr def cloud 
   do
      echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      echo "NOW PROCESSING ${var} FOR NEIGH RAD ${neigh_rad} METERS..."
      # NOTE: call Julia v1.4.1, now removed from `module avail`
      /global/home/groups/consultsw/sl-7.x86_64/modules/julia/1.4.1/bin/julia -p 32 /global/home/users/drewhart/seasonality/seasonal_asynchrony/src/asynch/calc/calc_asynch.jl --neigh_rad $neigh_rad --var $var > ch3_asynch_job_${var}_${neigh_rad}.jlout
      echo "FINISHED PROCESSING ${var} FOR NEIGH RAD ${neigh_rad} METERS."
   done
done
