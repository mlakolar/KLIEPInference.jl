#!/bin/bash

for m in "25" "50" "100"
do
  for sgn in "-1" "0" "1"
  do
    echo "sbatch ${m} ${sgn}..."
    sbatch --job-name=exp3_${m}_${sgn} --mem=4G sbatch_exp3 ${m} ${sgn}
    echo "... done"
  done
done
