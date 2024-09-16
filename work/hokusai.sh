#!/bin/bash
#
#SBATCH --job-name=test_mpi
#SBATCH --partition=mpc
#SBATCH --account=RB240033
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=5G
module load intel
export OMP_STACKSIZE=512M
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
srun --cpus-per-task=${SLURM_CPUS_PER_TASK} ./hormone
