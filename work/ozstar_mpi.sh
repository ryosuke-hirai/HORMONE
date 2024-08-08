#!/bin/bash
#
#SBATCH --job-name=test
#
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=32
#SBATCH --time=0:10:00
#SBATCH --mem-per-cpu=20G
module load intel-compilers/2023.0.0
module load gcccore/12.2.0
export OMP_STACKSIZE=512M
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
OMP_PLACES=cores
OMP_PROC_BIND=close
module load openmpi/4.1.5
srun ./hormone
