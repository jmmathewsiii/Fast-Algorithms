#!/bin/bash
#SBATCH --job-name=nbody-serial
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --time=0:05:00
#SBATCH --output=logs/serial-%j.out

module purge
module load gcc/11.2.0

cd $SLURM_SUBMIT_DIR
mkdir -p output scripts logs
./bin/serial
