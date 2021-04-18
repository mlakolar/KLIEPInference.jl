#!/bin/bash
SDIR="./res"
for m in "25" "50" "100"
do
  for sgn in "-1" "0" "1"
  do
    echo "batch inference ${m} ${sgn} ..."
    parallel --delay .2 julia experiment3.jl ${m} ${sgn} {1} $SDIR ::: {1..1000} 
    echo "... done"
  done
done
