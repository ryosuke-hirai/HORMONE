#!/bin/bash
#
#SBATCH --job-name=scaling_test
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=0:25:00
#SBATCH --mem=20G
module load intel-compilers/2023.0.0
module load gcccore/12.2.0
module load aocc/4.0.0
ulimit -s unlimited
export OMP_STACKSIZE=512M
rm scaling.dat
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=true
export OMP_DISPLAY_AFFINITY=true
export HORMONE_SCALING_TEST=true
for i in 64 32 16 8 4 2 1
do
 export OMP_NUM_THREADS=$i
 ./hormone
done
