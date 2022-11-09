#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH --cpus-per-task 8
#SBATCH --time 01:00:00

# Reference genome
RG=$1

# Species directory
SD=$2

# Working directory
WD=$3

# Sample directory
sample=$4

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.idxstats

exit 0