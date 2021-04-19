#!/bin/bash

#Change the scratch folder to somewhere you would like to store intermediate experimental results
SCRATCH_DIR=~/scratch

for graph in "chain1_25" "chain2_25"
do
  echo "sbatch $graph ..."
  sbatch --job-name=exp5_${graph} sbatch_exp5 ${graph} 15 150 300 $SCRATCH_DIR
  echo "... done"
done
for graph in "chain1_50" "chain2_50"
do
  echo "sbatch $graph ..."
  sbatch --job-name=exp5_${graph} sbatch_exp5 ${graph} 15 300 600 $SCRATCH_DIR
  echo "... done"
done
for graph in "tree1_25" "tree2_25"
do
  echo "sbatch $graph ..."
  sbatch --job-name=exp5_${graph} sbatch_exp5 ${graph} 2 150 300 $SCRATCH_DIR
  echo "... done"
done
for graph in "tree1_50" "tree2_50"
do
  echo "sbatch $graph ..."
  sbatch --job-name=exp5_${graph} sbatch_exp5 ${graph} 2 300 600 $SCRATCH_DIR
  echo "... done"
done
