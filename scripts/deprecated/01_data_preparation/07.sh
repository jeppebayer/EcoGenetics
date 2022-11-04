#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 8:00:00

# Mark duplicates and save output as .coordinatesort.bam. Also outputs some realted statistics
samtools markdup -@ 7 -s \
-d 100 \
-f people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.markdupstats \
-T people/Jeppe_Bayer/steps/temp/ \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.positionsort.bam \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.markdup.bam

exit 0