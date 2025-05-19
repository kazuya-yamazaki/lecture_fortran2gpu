#!/bin/bash
#PBS -q lecture-g
#PBS -l select=4:mpiprocs=8:ompthreads=9
#PBS -l walltime=00:10:00
#PBS -W group_list=gt00
#PBS -j oe

cd $PBS_O_WORKDIR

module purge
module load nvidia nv-hpcx                                                             

BIN=01_naive

mpirun -n ${MPI_PROC} ./${BIN} 


