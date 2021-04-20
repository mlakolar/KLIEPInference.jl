#!/bin/bash
SCRATCH_DIR="./res"

for m in "25" "50" "100"
do
  for sgn in "1"
  do
    for numChanges in "1" "3" "5"
    do
      for lbInd in {1..11}
      do
        echo "sbatch ${m} ${sgn} ${numChanges} ${lbInd} ..."
        # Change -j parameter to number of cores/threads in your CPU
        parallel --delay .2 -j 8 julia experiment4.jl ${m} ${sgn} ${numChanges} ${lbInd} {1} ${SCRATCH_DIR} ::: {1..1000}
        echo "... done"
      done
    done
  done
done