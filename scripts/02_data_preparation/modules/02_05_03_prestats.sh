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

# Creates coverage file for alignment
samtools coverage \
-o "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.coverage \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

exit 0