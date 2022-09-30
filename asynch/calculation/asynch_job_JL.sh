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

julia -p 32 /global/home/users/drewhart/seasonality/seasonal_asynchrony/asynch/calculation/calc_asynch.jl > ch3_asynch_job.jlout
