#!/bin/bash
# Job name:
#SBATCH --job-name=rasterr
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

module load python/3.7

# run the error-rasterization script
python /global/home/users/drewhart/seasonality/seasonal_asynchrony/asynch/src/anal/drivers/summ_results/rasterize_err.py > rast_err.pyout
