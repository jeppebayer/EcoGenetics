#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 8
#SBATCH --time 48:00:00

samtools view -b -@ 7 \
people/Jeppe_Bayer/steps/bwa/Orchesella_villosa?/FÅJ_C5/E100049659_L01_ANIcnqkR228559-607.sam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.bam