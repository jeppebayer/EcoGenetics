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

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ "$(("$cpus" - 1))" \
-F 3844 -q 20 \
-o "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ "$(("$cpus" - 1))" -b \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam.bai

exit 0