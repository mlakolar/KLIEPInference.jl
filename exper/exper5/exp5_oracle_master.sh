#!/bin/bash

expername=exp4_oracle

for graphtype in "chain_25_1"
do
  echo "sbatch $expername $graphtype ..."
  sbatch --job-name=${expername}_${graphtype} run_exp4_oracle.sh /home/bkim6/ising/ ${graphtype} 15 300
  echo "... done"
done

for graphtype in "chain_50_0" "chain_50_1"
do
  echo "sbatch $expername $graphtype ..."
  sbatch --job-name=${expername}_${graphtype} run_exp4_oracle.sh /home/bkim6/ising/ ${graphtype} 15 600
  echo "... done"
done

for graphtype in "tree_50_0" "tree_50_1"
do
  echo "sbatch $expername $graphtype ..."
  sbatch --job-name=${expername}_${graphtype} run_exp4_oracle.sh /home/bkim6/ising/ ${graphtype} 2 600
  echo "... done"
done

for graphtype in "chain_25_0"
do
  echo "sbatch $expername $graphtype ..."
  sbatch --job-name=${expername}_${graphtype} run_exp4_oracle.sh /home/bkim6/balanced/ ${graphtype} 15 300
  echo "... done"
done

for graphtype in "tree_25_0_2" "tree_25_1_1"
do
  echo "sbatch $expername $graphtype ..."
  sbatch --job-name=${expername}_${graphtype} run_exp4_oracle.sh /home/bkim6/balanced/ ${graphtype} 2 300
  echo "... done"
done
