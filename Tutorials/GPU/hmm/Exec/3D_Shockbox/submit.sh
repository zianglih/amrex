#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=run_mvc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=4000m
#SBATCH --time=01:00
#SBATCH --account=ramanvr0
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

module load cuda gcc openmpi

mpirun -np 2 HMM3d.gnu.TPROF.MPI.CUDA.ex inputs
