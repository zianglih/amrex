#!/bin/bash
#SBATCH --job-name=run_test
#SBATCH --time=01:00
#SBATCH --account=ramanvr0
#SBATCH --partition=gpu

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1g
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=closest
##SBATCH --gpu-bind=map_gpu:0,1
module load cuda
ncu --query-metrics