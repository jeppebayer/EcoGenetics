#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

for file in "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/*; do
    if [ ! -s "$file" ]; then
        rm -f "$file"
    fi
done