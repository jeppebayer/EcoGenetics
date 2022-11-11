#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

# Reference genome
RG=$1

# Species directory
SD=$2

# Working directory
WD=$3

# Sample directory
sample=$4

for file in "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/*; do
    if [ ! -s "$file" ]; then
        rm -f "$file"
    fi
done