#!/bin/bash


for set_i in "set1" "set2" "set3" "set4"
do
for delta in `seq -0.75 0.05 0.75`
do
    echo "sbatch  ${set_i} ${delta} ..."
    mkdir -p /scratch/mkolar/KLIEP/exp2_${set_i}_${delta}
    mkdir -p /scratch/mkolar/KLIEP/out/exp2_${set_i}_${delta}
    sbatch --job-name=exp2_${set_i}_${delta} sbatch_exp2 ${set_i} ${delta}
    echo "... done"
done
done
