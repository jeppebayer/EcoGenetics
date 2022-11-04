#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 6
#SBATCH --time 48:00:00

bwa mem -t 6 \
BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic \
BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_1.fq.gz BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_2.fq.gz \
> people/Jeppe_Bayer/steps/bwa/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sam