#!/bin/bash
#PBS -N Blastwave2D
#PBS -q openmp
#PBS -m n
#PBS -l select=1:ncpus=4
#PBS -l walltime=00:30:00
#PBS -j oe

module load intel/2024

#export OMP_NUM_THREADS=${PBS_NCPUS}
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cd $PBS_O_WORKDIR
./Simulation.x

