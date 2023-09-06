#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.idxstats

exit 0