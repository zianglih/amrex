#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=run_cns_sod
#SBATCH --time=05:00
#SBATCH --account=ramanvr0
#SBATCH --partition=gpu

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=2
#SBATCH --mem-per-cpu=8g
#SBATCH --gpu-bind=map_gpu:0,1

module purge

module load cuda cupti gcc openmpi

# For MPI:
# nsys profile --trace=cuda,nvtx,osrt,mpi --cuda-memory-usage=true mpirun -np 2 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs

# For a single CPU & GPU:
# nsys profile --trace=cuda,nvtx,osrt --cuda-memory-usage=true mpirun -np 1 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs

# For Nsight Compute:
ncu --nvtx --nvtx-include "main()" --target-processes all CNS3d.gnu.TPROF.MPI.CUDA.ex inputs amrex.fpe_trap_invalid=0
