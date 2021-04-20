#!/bin/bash
SCRATCH_DIR="./res"
for m in "25" "50" "100"
do
  for sgn in "-1" "0" "1"
  do
    echo "batch inference ${m} ${sgn} ..."
    # Change -j parameter to number of cores/threads in your CPU
    parallel --delay .2 -j 16 julia experiment3.jl ${m} ${sgn} {1} $SCRATCH_DIR ::: {1..1000} 
    echo "... done"
  done
done
