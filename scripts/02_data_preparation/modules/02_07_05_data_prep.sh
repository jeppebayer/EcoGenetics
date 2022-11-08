#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
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

# Creates coverage file for alignment
samtools stats -@ 7 \
"$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.stats

exit 0