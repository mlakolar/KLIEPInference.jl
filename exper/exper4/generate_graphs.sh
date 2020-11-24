#!/bin/bash

for numChanges in "1" "3" "5"
do
  for lbInd in {1..11}
  do
    echo "jl1 exp2_generate_graph.jl 100 1 ${numChanges} ${lbInd} ..."
    jl1 exp2_generate_graph.jl 100 1 ${numChanges} ${lbInd}
    echo "... done"
  done
done
