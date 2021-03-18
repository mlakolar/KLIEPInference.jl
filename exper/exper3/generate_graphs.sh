#!/bin/bash

for m in "25" "50" "100"
do
  for sgn in "-1" "0" "1"
  do
    echo "julia exp3_generate_graphs.jl ${m} ${sgn} ..."
    julia exp3_generate_graphs.jl ${m} ${sgn}
    echo "... done"
  done
done
