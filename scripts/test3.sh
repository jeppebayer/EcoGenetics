#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

RG=$1

path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')

srun echo "$path"

exit 0