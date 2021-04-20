#!/bin/bash
SCRATCH_DIR="/mnt/storage/home/sl9885/scratch/expr4"

for m in "25"
do
  for sgn in "1"
  do
    for numChanges in "1"
    do
      for lbInd in {1..11}
      do
        echo "sbatch ${m} ${sgn} ${numChanges} ${lbInd} ..."
        sbatch --job-name=exp4_${m}_${sgn}_${numChanges}_${lbInd} sbatch_exp4 ${m} ${sgn} ${numChanges} ${lbInd} $SCRATCH_DIR
        echo "... done"
      done
    done
  done
done
