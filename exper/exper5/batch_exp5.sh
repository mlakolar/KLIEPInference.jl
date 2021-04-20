#!/bin/bash

#Change the scratch folder to somewhere you would like to store intermediate experimental results
SCRATCH_DIR=~/scratch

for graph in "chain1_25" "chain2_25"
do
  echo "$graph ..."
  parallel --delay .2 -j 16 julia experiment5.jl ${graph} 15 150 300 {1} $SCRATCH_DIR ::: {1..1000} 
  echo "... done"
done
for graph in "chain1_50" "chain2_50"
do
  echo "$graph ..."
  parallel --delay .2 -j 16 julia experiment5.jl ${graph} 15 300 600 {1} $SCRATCH_DIR ::: {1..1000} 
  echo "... done"
done
for graph in "tree1_25" "tree2_25"
do
  echo "$graph ..."
  parallel --delay .2 -j 16 julia experiment5.jl ${graph} 2 150 300 {1} $SCRATCH_DIR ::: {1..1000} 
  echo "... done"
done
for graph in "tree1_50" "tree2_50"
do
  echo "$graph ..."
  parallel --delay .2 -j 16 julia experiment5.jl ${graph} 2 300 600 {1} $SCRATCH_DIR ::: {1..1000} 
  echo "... done"
done
