#!/bin/bash
# Job name:
#SBATCH --job-name=EOF_clust
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

echo "CALCULATING EOFs AND CLUSTERS FOR GLOBAL NIRv LSP MAP..."
python /global/home/users/drewhart/seasonality/seasonal_asynchrony/src/phen/anal/div/calc_global_LSP_EOFs_and_clusts.py > LSP_EOFs.pyout
echo "FINISHED CALCULATING EOFS AND CLUSTERS."
