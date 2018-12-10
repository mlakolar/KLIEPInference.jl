#!/bin/bash

for numChanges in "1" "3" "5"
do
  for lbInd in {1..11}
  do
    echo "sbatch 100 1 ${numChanges} ${lbInd} ..."
    sbatch --job-name=exp2_100_1_${numChanges}_${lbInd}_500_500 sbatch_exp2 100 1 ${numChanges} ${lbInd} 500 500
    echo "... done"
  done
done
