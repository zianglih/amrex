#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=run_cns_sod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2048m
#SBATCH --time=05:00
#SBATCH --account=ramanvr0
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

module load cuda cupti gcc openmpi

nsys profile --trace=cuda,nvtx,osrt --backtrace=dwarf --cuda-memory-usage=true mpirun -np 1 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs
# nvprof mpirun -np 1 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs
# nsys profile -c nvtx -p "main@*" -e NSYS_NVTX_PROFILER_REGISTER_ONLY=0 --trace=cuda,nvtx,osrt --backtrace=dwarf --gpu-metrics-device=all --gpu-metrics-frequency=20000 mpirun -np 1 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs
# mpirun -np 1 ncu -o kernels --target-processes all -c 2 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs amrex.fpe_trap_invalid=0
# ncu -k 5 --target-processes all CNS3d.gnu.TPROF.MPI.CUDA.ex inputs amrex.fpe_trap_invalid=0
