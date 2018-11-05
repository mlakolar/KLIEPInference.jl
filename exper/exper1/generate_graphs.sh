#!/bin/bash

for p in "25" "50" "100"
do
  for sgn in "-1" "0" "1"
  do
    echo "jl1 exp1_generate_graph.jl ${p} ${sgn} ..."
    jl1 exp1_generate_graph.jl ${p} ${sgn}
    echo "... done"
  done
done
