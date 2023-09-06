#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ 7 \
-F 0x404 -f 0x2 -q 20 \
-o people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.markdup.bam