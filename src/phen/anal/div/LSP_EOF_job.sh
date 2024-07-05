#!/bin/bash
# Job name:
#SBATCH --job-name=LSP_EOF
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3_bigmem
#
# Wall clock limit:
#SBATCH --time=12:00:00
#
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Commands to run:

stdbuf -i0 -o0 -e0 command

module load python/3.7

echo "CALCULATING EOFs FOR NIRv-BASED LSP MAP..."
python /global/home/users/drewhart/seasonality/seasonal_asynchrony/phen/analysis/div/calc_LSP_EOFs.py > LSP_EOFs.pyout
echo "FINISHED CALCULATING EOFS."
