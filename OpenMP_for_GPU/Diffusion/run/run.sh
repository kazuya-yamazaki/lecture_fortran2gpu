#!/bin/bash
#PBS -q lecture-g
#PBS -l select=4:mpiprocs=1
#PBS -l walltime=00:03:00
#PBS -W group_list=gt00
#PBS -j oe

cd $PBS_O_WORKDIR

module purge
module load nvidia nv-hpcx

mpirun -n 4 ./wrapper.sh ./ompgpu_unified
mpirun -n 4 ./wrapper.sh ./ompgpu
