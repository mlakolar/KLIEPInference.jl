#!/bin/bash

for m in "25"
do
  for sgn in "0"
  do
    for numChanges in "1" "3" "5"
    do
      for lbInd in {1..11}
      do
        echo "sbatch ${m} ${sgn} ${numChanges} ${lbInd} ..."
        sbatch --job-name=exp4_${m}_${sgn}_${numChanges}_${lbInd} sbatch_exp4 ${m} ${sgn} ${numChanges} ${lbInd}
        echo "... done"
      done
    done
  done
done
