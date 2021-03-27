#!/bin/bash

for m in "25"
do
  for sgn in "0"
  do
    for numChanges in "1" "3" "5"
    do
      for lbInd in {1..11}
      do
        echo "julia exp4_generate_graph.jl ${m} ${sgn} ${numChanges} ${lbInd} ..."
        julia exp4_generate_graph.jl ${m} ${sgn} ${numChanges} ${lbInd}
        echo "... done"
      done
    done
  done
done
