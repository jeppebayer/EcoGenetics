#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Align sample to reference genome and pipe output
bwa mem -t 8 \
-o people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.raw.sam \
BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic \
BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_1.fq.gz BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_2.fq.gz

exit 0