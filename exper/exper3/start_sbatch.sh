#!/bin/bash

## This script will start 9 nodes, for 9 differernt settings of m and sgn.
ScratchDIR="/mnt/storage/home/sl9885/scratch/expr3"

for m in "25"
do
  for sgn in "-1"
  do
    echo "sbatch ${m} ${sgn} ..."
    sbatch --job-name=exp3_${m}_${sgn} sbatch_exp3 ${m} ${sgn} $ScratchDIR
    echo "... done"
  done
done
