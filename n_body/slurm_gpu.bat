#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=nbody-cuda
#SBATCH --partition=atesting_a100
#SBATCH --gres=gpu:1
#SBATCH --qos=testing
#SBATCH --ntasks=8
#SBATCH --time=00:05:00
#SBATCH --output=logs/cuda-%j.out

module purge
module load cuda gcc/11.2.0

cd $SLURM_SUBMIT_DIR
mkdir -p output scripts logs
./bin/cuda
