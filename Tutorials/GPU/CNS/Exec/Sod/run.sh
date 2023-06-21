#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=run_cns_sod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=8000m
#SBATCH --time=05:00
#SBATCH --account=ramanvr0
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

module load cuda gcc openmpi

nsys profile mpirun -np 4 CNS3d.gnu.TPROF.MPI.CUDA.ex inputs
