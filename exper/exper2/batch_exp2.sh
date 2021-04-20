#!/bin/bash

#Change the scratch folder to somewhere you would like to store intermediate experimental results
SCRATCH_DIR=~/scratch

for set_i in "set1" "set2" "set3" "set4"
do
for delta in `seq -0.75 0.05 0.75`
do
    echo "${set_i} ${delta} ..."
    mkdir -p ${SCRATCH_DIR}/exp2_${set_i}_${delta}
    mkdir -p ${SCRATCH_DIR}/out/exp2_${set_i}_${delta}
    # Change -j parameter to number of cores/threads in your CPU
    parallel --delay .2 -j 16 julia experiment2.jl ${set_i} ${delta} {1} ${SCRATCH_DIR} ::: {1..10}
    echo "... done"
done
done