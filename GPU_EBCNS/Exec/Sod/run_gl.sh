#!/bin/bash

#SBATCH --job-name=run_gpu_ebcns_sod
#SBATCH --time=10:00
#SBATCH --account=ramanvr0
#SBATCH --partition=gpu

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=16g
#SBATCH --gpus-per-node=2
#SBATCH --gpu-bind=closest
##SBATCH --gpu-bind=map_gpu:0,1

module load cuda gcc openmpi

# For Nsight Systems:
nsys profile --trace=cuda,nvtx,osrt,mpi --cuda-memory-usage=true mpirun -np 2 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs

# For Nsight Compute:

# Currently, no kernels will be profiled from this command
# ncu --nvtx --nvtx-include "CNS()" --target-processes all CNS3d.gnu.TPROF.MPI.CUDA.ex inputs amrex.fpe_trap_invalid=0

# Get this command from right clicking in Nsight Systems
# ncu --target-processes all --nvtx --call-stack -o compute_dSdT_kernels_metrics --kernel-name launch_global --launch-skip 446 --launch-count 31 mpirun -np 2 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs

