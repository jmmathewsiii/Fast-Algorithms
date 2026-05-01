#!/bin/bash
#SBATCH --job-name=nbody-serial
#SBATCH --partition=atesting
#SBATCH --ntasks=1
#SBATCH --qos=testing
#SBATCH --cpus-per-task=1
#SBATCH --time=0:05:00
#SBATCH --output=logs/serial-%j.out

module purge
module load gcc/11.2.0

cd $SLURM_SUBMIT_DIR
mkdir -p output output/lyapunov scripts scripts/lyapunov logs

# Edit N / n_iter as needed. Mode is "speed" so timings are printed.
N=1000
NITER=1000

# Run direct first so BH has a baseline to validate against.
./bin/serial_direct <<EOF
$N
$NITER
speed
EOF

./bin/serial_bh <<EOF
$N
$NITER
speed
EOF
