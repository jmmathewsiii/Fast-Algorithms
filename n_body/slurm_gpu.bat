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
mkdir -p output output/lyapunov scripts scripts/lyapunov logs

# Edit N / n_iter as needed. Mode is "speed" so timings are printed.
N=10000
NITER=4000

./bin/cuda_direct <<EOF
$N
$NITER
speed
EOF
