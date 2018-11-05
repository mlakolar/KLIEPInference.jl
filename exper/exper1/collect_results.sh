#!/bin/bash

echo "spKLIEP"
jl1 exp1_cc.jl 25 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1
echo "oracle KLIEP"
jl1 exp1_cc.jl 25 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1/oracle

echo "spKLIEP"
jl1 exp1_cc.jl 50 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1
echo "oracle KLIEP"
jl1 exp1_cc.jl 50 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1/oracle
