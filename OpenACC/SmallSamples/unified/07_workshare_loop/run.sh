#!/bin/bash
#PBS -q lecture-mig
#PBS -l select=1:mpiprocs=1
#PBS -l walltime=00:01:00
#PBS -W group_list=gt00
#PBS -j oe

cd $PBS_O_WORKDIR

module load nvidia

./run
