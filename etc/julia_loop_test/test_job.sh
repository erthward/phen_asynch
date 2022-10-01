#!/bin/bash
# Job name:
#SBATCH --job-name=testloop
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=1:00:00
#
## Commands to run:

module load julia

rm /global/home/users/drewhart/testloop.jlout
for n in 2 7 1 5 23 3
do
	start=`date +%s`

	julia -p 4 /global/home/users/drewhart/test_script.jl -x $n  >> testloop.jlout
	end=`date +%s`
	runtime=$((end-start))
	echo "DONE WITH ${n}!"
	echo "(took ${runtime} s)"
done
