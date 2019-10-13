#!/bin/bash

for expername in "exp4_x" "exp4_y" "exp4_sym"
do
  for graphtype in "chain_25_0" "chain_25_1"
  do
    echo "sbatch $expername $graphtype ..."
    sbatch --job-name=${expername}_${graphtype} run_exp4.sh ${expername} ${graphtype} 15
    echo "... done"
  done
  for graphtype in "chain_50_0" "chain_50_1"
  do
    echo "sbatch $expername $graphtype ..."
    sbatch --job-name=${expername}_${graphtype} run_exp4.sh ${expername} ${graphtype} 15
    echo "... done"
  done
  for graphtype in "tree_25_0_2" "tree_25_1_1"
  do
    echo "sbatch $expername $graphtype ..."
    sbatch --job-name=${expername}_${graphtype} run_exp4.sh ${expername} ${graphtype} 2
    echo "... done"
  done
  for graphtype in "tree_50_0" "tree_50_1"
  do
    echo "sbatch $expername $graphtype ..."
    sbatch --job-name=${expername}_${graphtype} run_exp4.sh ${expername} ${graphtype} 2
    echo "... done"
  done
done
