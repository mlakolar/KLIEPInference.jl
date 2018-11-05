#!/bin/bash

for p in "25" "50" "100"
do
#  for sgn in "-1" "0" "1"
#  do
    echo "sbatch ${p} ..."
    sbatch --job-name=exp1_${p}_1_500_500 sbatch_exp1 ${p} 1 500 500
    echo "... done"
#  done
done
