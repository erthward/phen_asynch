#!/bin/bash
# Job name:
#SBATCH --job-name=window
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
## Command(s) to run (example):

stdbuf -i0 -o0 -e0 command

module load r r-packages r-spatial

Rscript --vanilla ./calc_circular_moving_window_chelsa_rasters.r > moving_window.Rout
