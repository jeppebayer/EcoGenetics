#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 24:00:00

# Reference genome
RG=$1

# Species directory
SD=$2

# Working directory
WD=$3

# Sample directory
sample=$4

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ 7 \
-F 3844 -q 20 \
-o "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
"$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ 7 -b \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam.bai

# File removal
rm -f \
"$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \

exit 0