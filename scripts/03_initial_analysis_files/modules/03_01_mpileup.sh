#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
dataprep=$6 # Path to script location
algo=$7 # Chosen algorithm

# Make mpileup file
samtools mpileup \
-C 0 \
-o "$sample"/"$(basename "$sample")".pileup \
-f "$RG" \
"$sample"/"$(basename "$sample")"_filtered.bam

exit 0