#!/bin/bash
#
#SBATCH --job-name=test_mpi
#
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=32
#SBATCH --time=0:10:00
#SBATCH --mem-per-cpu=10G
module load gcc/12.3.0
module load openmpi/4.1.5
export OMP_STACKSIZE=512M
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
OMP_PLACES=cores
OMP_PROC_BIND=close
srun ./hormone
