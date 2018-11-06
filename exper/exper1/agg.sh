#!/bin/bash


for p in "25" "50" "100"
do
  for sgn in "-1" "0" "1"
  do
    echo "sbatch ${p} ${sgn} ..."
    jl1 aggregate_results.jl ${p} ${sgn} 500 500 /scratch/midway2/mkolar/KLIEP/exp1/oracle agg_oracle
    echo "... done"
  done
done



