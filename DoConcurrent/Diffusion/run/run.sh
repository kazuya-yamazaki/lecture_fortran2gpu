#!/bin/bash
#PBS -q lecture-g
#PBS -l select=4:mpiprocs=1
#PBS -l walltime=00:02:00
#PBS -W group_list=gt00
#PBS -j oe

cd $PBS_O_WORKDIR

module purge
module load nvidia nv-hpcx

mpirun -n ${MPI_PROC} ./wrapper.sh ./stdpar
#mpirun -n ${MPI_PROC} nsys profile -f true -t cuda,nvtx,openacc,mpi -o report_%q{OMPI_COMM_WORLD_RANK} ./wrapper.sh ./stdpar
