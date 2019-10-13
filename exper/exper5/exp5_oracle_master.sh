#!/bin/bash

for graphtype in "chain_25_0" "chain_25_1" "chain_50_0" "chain_50_1"
do
  echo "sbatch exp5_oracle $graphtype ..."
  sbatch --job-name=exp5_oracle_${graphtype} run_exp5_oracle.sh ${graphtype} 15
  echo "... done"
done
for graphtype in "tree_25_0_2" "tree_25_1_1" "tree_50_0" "tree_50_1"
do
  echo "sbatch exp5_oracle $graphtype ..."
  sbatch --job-name=exp5_oracle_${graphtype} run_exp5_oracle.sh ${graphtype} 2
  echo "... done"
done
