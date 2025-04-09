#!/bin/bash
#
#SBATCH --job-name=hormone_lc
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=50:10:00
#SBATCH --mem=5G
#module load intel-compilers/2023.0.0
#module load gcccore/12.2.0
export OMP_STACKSIZE=512M
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
ulimit -s unlimited
OMP_PLACES=cores
OMP_PROC_BIND=close
./lightcurve
