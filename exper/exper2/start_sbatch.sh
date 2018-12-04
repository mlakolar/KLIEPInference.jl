#!/bin/bash



for p in "25" "50" "100"
do
  for numChanges in "1" "3" "5"
  do
    for lbInd in {1..11}
    do
      echo "sbatch ${p} 1 ${numChanges} ${lbInd} ..."
      sbatch --job-name=exp2_${p}_1_${numChanges}_${lbInd}_500_500 sbatch_exp2 ${p} 1 ${numChanges} ${lbInd} 500 500
      echo "... done"
    done
  done
done
