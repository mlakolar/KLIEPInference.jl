#!/bin/bash



for p in "25" "50" 
do
  for numChanges in "4" "8" 
  do
    for lbInd in "1" "2" "3" 
    do
      echo "sbatch ${p} ${numChanges} ${lbInd} ..."
      jl1 aggregate_results.jl ${p} 1 ${numChanges} ${lbInd} 500 500 /scratch/midway2/mkolar/KLIEP/exp2 agg_empBoot
      echo "... done"
    done
  done
done



