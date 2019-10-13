#!/bin/bash

for expername in "x" "y"
do
  for graphtype in "chain_25_0" "chain_25_1" "chain_50_0" "chain_50_1"
  do
    echo "sbatch sanity check $expername $graphtype ..."
    sbatch --job-name=sanitycheck_${expername}_${graphtype} run_sanitycheck.sh ${expername} ${graphtype} 15
    echo "... done"
  done
  for graphtype in "tree_25_0_2" "tree_25_1_1" "tree_50_0" "tree_50_1"
  do
    echo "sbatch sanity check $expername $graphtype ..."
    sbatch --job-name=sanitycheck_${expername}_${graphtype} run_sanitycheck.sh ${expername} ${graphtype} 2
    echo "... done"
  done
done
