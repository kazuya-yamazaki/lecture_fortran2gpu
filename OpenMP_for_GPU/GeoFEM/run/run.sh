#!/bin/bash
#PBS -q lecture-g
#PBS -l select=4:mpiprocs=1
#PBS -l walltime=00:02:00
#PBS -W group_list=gt00
#PBS -j oe

cd $PBS_O_WORKDIR

for exe in Omp_U_a Omp_U_c
do
    echo ${exe}
    mpirun -n ${MPI_PROC} ./wrapper.sh ./${exe}
done

