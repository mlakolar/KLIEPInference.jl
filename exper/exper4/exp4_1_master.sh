#!/bin/bash

for expername in "exp4_y" "exp4_x" "exp4_sym"
do
  for graphtype in "chain_25_1"
  do
    echo "sbatch $expername $graphtype ..."
    sbatch --job-name=${expername}_${graphtype} run_exp4.sh ${expername} /home/bkim6/ising/ ${graphtype} 15 300
    echo "... done"
  done
  for graphtype in "chain_50_0" "chain_50_1"
  do
    echo "sbatch $expername $graphtype ..."
    sbatch --job-name=${expername}_${graphtype} run_exp4.sh ${expername} /home/bkim6/ising/ ${graphtype} 15 600
    echo "... done"
  done
  for graphtype in "tree_50_0" "tree_50_1"
  do
    echo "sbatch $expername $graphtype ..."
    sbatch --job-name=${expername}_${graphtype} run_exp4.sh ${expername} /home/bkim6/ising/ ${graphtype} 2 600
    echo "... done"
  done
done
