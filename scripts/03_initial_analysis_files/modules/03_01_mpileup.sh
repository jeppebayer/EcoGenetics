#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

# Reference genome
RG=$1

# Species directory
SD=$2

# Working directory
WD=$3

# Sample directory
sample=$4

# Make mpileup file
samtools mpileup \
-C 0 \
-o "$sample"/"$(basename "$sample")".mpileup \
-f "$RG" \
"$sample"/"$(basename "$sample")"_filtered.bam

exit 0