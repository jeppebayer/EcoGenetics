#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates stat file for alignment
samtools stats -@ "$((cpus - 1))" \
-c 1,1000,1 \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.stats

exit 0