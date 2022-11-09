#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Reference genome
RG=$1

# Species directory
SD=$2

# Working directory
WD=$3

# Sample directory
sample=$4

# Creates stat file for alignment
samtools stats -@ 7 \
-c 1,1000,1 \
"$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.stats

exit 0