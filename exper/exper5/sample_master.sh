#!/bin/bash

for graphtype in "chain_25_0" "chain_25_1" "tree_25_0_2" "tree_25_1_1"
do
  echo "sbatch sample $graphtype ..."
  sbatch --job-name=sample_${graphtype} sample_ising.sh ${graphtype} 300
  echo "... done"
done
for graphtype in "chain_50_0" "chain_50_1" "tree_50_0" "tree_50_1"
do
  echo "sbatch sample $graphtype ..."
  sbatch --job-name=sample_${graphtype} sample_ising.sh ${graphtype} 600
  echo "... done"
done
