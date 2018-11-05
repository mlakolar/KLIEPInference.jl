#!/bin/bash

for p in "25" "50" "100"
do
  for numChanges in "4" "8" 
  do
    for lbInd in "1" "2" "3" 
    do
      echo "jl1 exp2_generate_graph.jl ${p} 1 ${numChanges} ${lbInd} ..."
      jl1 exp2_generate_graph.jl ${p} 1 ${numChanges} ${lbInd} 
      echo "... done"
    done
  done
done
