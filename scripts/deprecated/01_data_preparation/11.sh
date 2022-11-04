#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 4:00:00

# Creates coverage file for alignment
samtools coverage \
-o people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.coverage \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.markdup.bam

exit 0