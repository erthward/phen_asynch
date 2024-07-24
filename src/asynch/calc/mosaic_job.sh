#!/bin/bash
# Job name:
#SBATCH --job-name=mosaic
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=5:00:00
#
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Commands to run:

stdbuf -i0 -o0 -e0 command

module load python/3.7

python /global/home/users/drewhart/seasonality/seasonal_asynchrony/asynch/src/calc/mosaic_all_results.py > mosaic_all_results.pyout
