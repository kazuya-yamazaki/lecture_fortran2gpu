#!/bin/bash
#PBS -q lecture-g
#PBS -l select=4:mpiprocs=1
#PBS -l walltime=00:02:00
#PBS -W group_list=gt00
#PBS -j oe

cd $PBS_O_WORKDIR

module purge
module load nvidia nvmpi                                                             

#export NVIDIA_ACC_TIME=1

BIN=00_work

mpirun -n ${MPI_PROC} ./wrapper.sh ./${BIN} 


