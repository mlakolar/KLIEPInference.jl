#!/bin/bash

## m = 25
for sgn in "-1" "0" "1"
do
  echo "sbatch 25 ${sgn}..."
  sbatch --job-name=exp3_25_${sgn} sbatch_exp3 25 ${sgn}
  echo "... done"
done

## m = 50
for sgn in "-1" "0" "1"
do
  echo "sbatch 50 ${sgn}..."
  sbatch --job-name=exp3_50_${sgn} sbatch_exp3 50 ${sgn}
  echo "... done"
done

## m = 100
for sgn in "-1" "0" "1"
do
  echo "sbatch 100 ${sgn}..."
  sbatch --job-name=exp3_100_${sgn} --mem=4G sbatch_exp3 100 ${sgn}
  echo "... done"
done
