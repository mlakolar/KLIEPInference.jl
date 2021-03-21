#!/bin/bash

## call this file as
##   sbatch --job-name=exp5_graph sbatch_exp5 graph idx nx ny

#SBATCH --partition=broadwl
#SBATCH --time=1-12:00:00
#SBATCH --ntasks=168

module load parallel
module load hdf5

ERR="--error=/scratch/midway2/byolkim/$SLURM_JOB_NAME/$SLURM_JOB_NAME.err.{1}"
OUT="--output=/scratch/midway2/byolkim/$SLURM_JOB_NAME/$SLURM_JOB_NAME.out.{1}"

mkdir -p /scratch/midway2/byolkim/$SLURM_JOB_NAME/

parallel --delay .2 -j $SLURM_NTASKS --joblog __$SLURM_JOB_NAME.runtask.log --resume srun --exclusive -N1 -n1 $ERR $OUT /project2/mkolar/julia-1.5.0/bin/julia experiment5.jl ${2} ${3} ${4} ${5} {1} ::: {1..1000}
