#!/bin/bash
#PBS -N Blastwave2D
#PBS -q openmp
#PBS -m n
#PBS -l ncpu=4
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID}

module load intel/2024

export OMP_NUM_THREADS=$PBS_NCPUS
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cd $PBS_O_WORKDIR
./Simulation.x


