#!/bin/bash
#PBS -N Blastwave2D
#PBS -q small
#PBS -m n
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
cd $PBS_O_WORKDIR
./Simulation.x


