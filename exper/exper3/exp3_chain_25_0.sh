#!/bin/bash

#SBATCH --account=pi-mkolar

#SBATCH --job-name=single_chain_25_0

#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --time=2-12:00:00
#SBATCH --exclusive

module load julia/1.1.1

mkdir -p /scratch/bkim6/single/chain_25_0

ERR="--error=/scratch/bkim6/single/chain_25_0/$SLURM_JOB_ID.err.{1}"
OUT="--output=/scratch/bkim6/single/chain_25_0/$SLURM_JOB_ID.out.{1}"

parallel --delay .2 -j $SLURM_NTASKS --joblog __$SLURM_JOB_NAME.runtask.log --resume srun --exclusive -N1 -n1 $ERR $OUT julia experiment3.jl chain_25_0 15 {1} 0.316 150 300 ::: {1..1000}
