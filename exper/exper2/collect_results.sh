#!/bin/bash

for p in "25" "50"
do
  for numChanges in "1" "3" "5"
  do
    for lbInd in {1..11}
    do
      echo "spKLIEP ${p} 1 ${numChanges} ${lbInd} ..."
      /project/mkolar/julia-1.0.2/bin/julia exp2_cc.jl ${p} 1 ${numChanges} ${lbInd} 500 500 /scratch/midway2/byolkim/power ${1}
    done
  done
done
