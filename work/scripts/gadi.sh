#!/bin/bash

#PBS -P di94
#PBS -q normal
#PBS -l ncpus=768
#PBS -l mem=400gb
#PBS -l walltime=24:00:00
#PBS -l wd
#module purge
#module load pbs
#module load intel-compiler-llvm/2024.2.1
module load intel-mpi
#module load openmpi
export OMP_STACKSIZE=512M
export OMP_NUM_THREADS=1
OMP_PLACES=cores
OMP_PROC_BIND=close
mpiprocs=$(( $PBS_NCPUS / $OMP_NUM_THREADS ))
echo ${mpiprocs}
# Run the MPI+OpenMP hybrid job
mpirun -np ${mpiprocs} ./hormone
