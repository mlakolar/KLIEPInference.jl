#!/bin/bash

for p in "25" "50" "100"
do
  for sgn in "-1" "0" "1"
  do
    echo "sbatch ${p} ${sgn} ..."
    sbatch --job-name=oracle_exp1_${p}_${sgn}_500_500 sbatch_oracle_exp1 ${p} ${sgn} 500 500
    echo "... done"
  done
done
