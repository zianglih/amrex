#!/bin/bash

#SBATCH --job-name=run_ebcns_sod
#SBATCH --time=30:00
#SBATCH --account=ramanvr0
#SBATCH --partition=gpu

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=32g
#SBATCH --gpus-per-node=2
##SBATCH --gpu-bind=per_task:1
##SBATCH --gpu-bind=closest
#SBATCH --gpu-bind=map_gpu:0,1

module load cuda gcc openmpi

# For MPI:
nsys profile --trace=cuda,nvtx,osrt,mpi --cuda-memory-usage=true mpirun -np 2 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs

# For a single CPU & GPU:
# nsys profile --trace=cuda,nvtx,osrt --cuda-memory-usage=true mpirun -np 1 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs

# For Nsight Compute:
# Currently, no kernels will be profiled
# ncu --nvtx --nvtx-include "CNS()" --target-processes all CNS3d.gnu.TPROF.MPI.CUDA.ex inputs amrex.fpe_trap_invalid=0
