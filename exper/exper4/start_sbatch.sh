#!/bin/bash
ScratchDIR="/mnt/storage/home/sl9885/scratch/expr4"

for m in "25" "50" "100"
do
  for sgn in "1"
  do
    for numChanges in "1" "3" "5"
    do
      for lbInd in {1..11}
      do
        echo "sbatch ${m} ${sgn} ${numChanges} ${lbInd} ..."
        sbatch --job-name=exp4_${m}_${sgn}_${numChanges}_${lbInd} sbatch_exp4 ${m} ${sgn} ${numChanges} ${lbInd} $ScratchDIR
        echo "... done"
      done
    done
  done
done
