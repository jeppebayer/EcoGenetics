#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 8
#SBATCH --time 48:00:00

samtools sort -@ 7 \
people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.sorted.bam