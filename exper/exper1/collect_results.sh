#!/bin/bash

for p in "25" "50" "100"
do
    echo "============================================================="
    echo ""
    echo "spKliep ${1} ::: p = ${p}"
    jl1 exp1_cc.jl ${p} 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1 ${1}
    echo ""
    echo ""
    echo "oracleKliep ${1}"
    jl1 exp1_cc.jl ${p} 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1/oracle ${1}
    echo ""
    echo "============================================================="
done

