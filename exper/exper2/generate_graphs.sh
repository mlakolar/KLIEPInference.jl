#!/bin/bash

for p in "25" "50" "100"
do
  for numChanges in "1" "3" "5"
  do
    for lbInd in {1..11}
    do
      echo "jl102 exp2_generate_graph.jl ${p} 1 ${numChanges} ${lbInd} ..."
      /project/mkolar/julia-1.0.2/bin/julia exp2_generate_graph.jl ${p} 1 ${numChanges} ${lbInd}
      echo "... done"
    done
  done
done
