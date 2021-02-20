#!/bin/bash

for graphtype in "tree_25_0_2" "tree_25_1_1"
do
  echo "sbatch sample $graphtype ..."
  sbatch --job-name=sample_${graphtype} sample_ising.sh ${graphtype} 300
  echo "... done"
done
