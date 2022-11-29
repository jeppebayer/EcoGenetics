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

# Splits file into n parts without splitting lines. Each part gets samplename + number (00-99) + extension (.part)
split -n l/"$n" -d \
--additional-suffix .pileup \
"$sample"/"$(basename "$sample")".pileup \
"$WD"/temp/"$(basename "$sample")"_part

exit 0