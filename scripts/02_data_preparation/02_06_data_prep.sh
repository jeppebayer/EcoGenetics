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
-F 0x404 -f 0x2 -q 20 \
-o filtered.bam \
markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ 7 -b \
filtered.bam \
> filtered.bam.bai

exit 0