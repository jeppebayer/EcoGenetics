#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Creates bai index for alignment
samtools index -@ 7 -b \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.markdup.bam.bai

exit 0