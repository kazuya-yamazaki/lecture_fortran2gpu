#!/bin/bash                                                                                                               

gr=${OMPI_COMM_WORLD_RANK:-0}
lr=${OMPI_COMM_WORLD_LOCAL_RANK:-0}

NUM_GPUS=${NUM_GPUS:-8}

gpuid=`expr ${lr} \% ${NUM_GPUS}`
export CUDA_VISIBLE_DEVICES=${gpuid}
#export PGI_ACC_DEVICE_NUM=${gpuid}                                                                                       

# Launch kernels with "omp nowait" to the same stream                                                                     
export NV_OMP_AUTO_STREAMS=FALSE

nodeid=${gpuid}

if [ ${lr} == -1 ]; then
    nsys profile $@
    # gdb $@                                                                                                              
else
    numactl -N ${nodeid} $@
fi
