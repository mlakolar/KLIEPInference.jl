#!/bin/bash

echo "spKLIEP ${1}"
/project/mkolar/julia-1.0.2/bin/julia exp1_cc.jl 25 1 500 500 /scratch/midway2/byolkim/power ${1}
#echo "oracle KLIEP ${1}"
#jl1 exp1_cc.jl 25 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1/oracle ${1}

echo "spKLIEP ${1}"
/project/mkolar/julia-1.0.2/bin/julia exp1_cc.jl 50 1 500 500 /scratch/midway2/byolkim/power ${1}
#echo "oracle KLIEP ${1}"
#jl1 exp1_cc.jl 50 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1/oracle ${1}

#echo "spKLIEP"
#jl1 exp1_cc.jl 100 1 500 500 /scratch/midway2/byolkim/power
# echo "oracle KLIEP"
# jl1 exp1_cc.jl 100 1 500 500 /scratch/midway2/mkolar/KLIEP/exp1/oracle
