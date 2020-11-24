#!/bin/bash

echo "sbatch sanity check chain_25_0 ..."
sbatch --job-name=sanitycheck_chain_25_0 run_sanitycheck.sh chain_25_0 15
echo "... done"
