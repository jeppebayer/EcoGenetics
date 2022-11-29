#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
data=$6 # Path to data location in WD
n=$7 # Number of parts
script_path=$8 # Path to script location

# Make pileup file
samtools mpileup \
-C 0 \
-o "$sample"/"$(basename "$sample")".pileup \
-f "$RG" \
"$sample"/"$(basename "$sample")"_filtered.bam

exit 0