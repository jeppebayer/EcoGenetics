#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 8
#SBATCH --time 30:00:00

samtools view -b -@ 7 \
-F 0x4 -f 0x2 -q 20 \
people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/filtered/E100049659_L01_ANIcnqkR228559-607.bam