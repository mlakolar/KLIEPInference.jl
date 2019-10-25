#!/bin/bash

#SBATCH --account=pi-mkolar

## Call this with sbatch --job-name=${expername}_${graphtype} run_exp4.sh ${expername} ${graphdir} ${graphtype} ${idx} ${nx}

#SBATCH --partition=gpu
#SBATCH --ntasks=28
#SBATCH --time=48:00:00
#SBATCH --exclusive

ERR="--error=/scratch/bkim6/$SLURM_JOB_NAME/$SLURM_JOB_NAME.err.{1}"
OUT="--output=/scratch/bkim6/$SLURM_JOB_NAME/$SLURM_JOB_NAME.out.{1}"

mkdir -p /scratch/bkim6/$SLURM_JOB_NAME/

module load julia/1.1.1

parallel --delay .2 -j $SLURM_NTASKS --joblog __$SLURM_JOB_NAME.runtask.log --resume srun --exclusive -N1 -n1 $ERR $OUT julia ${1}.jl ${2} ${3} ${4} ${5} ${5} {1} ::: {1..1000}
