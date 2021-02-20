#!/bin/bash

for m in "25" "50" "100"
do
  for sgn in "-1" "0" "1"
  do
    echo "jl exp3_generate_graph.jl ${m} ${sgn} ..."
    jl exp3_generate_graph.jl ${m} ${sgn}
    echo "... done"
  done
done
