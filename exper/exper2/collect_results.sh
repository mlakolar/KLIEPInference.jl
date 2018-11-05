#!/bin/bash




for p in "25" "50" "100"
do
  for numChanges in "4" "8" 
  do
    for lbInd in "1" "2" "3" 
    do
      echo "============================================================="
      echo ""
      echo "spKliep ${1}"
      jl1 exp2_cc.jl ${p} 1 ${numChanges} ${lbInd} 500 500 /scratch/midway2/mkolar/KLIEP/exp2 ${1}
      echo ""
      echo ""
      echo "oracleKliep ${1}"
      jl1 exp2_cc.jl ${p} 1 ${numChanges} ${lbInd} 500 500 /scratch/midway2/mkolar/KLIEP/exp2/oracle ${1}
      echo ""
      echo "============================================================="
    done
  done
done


