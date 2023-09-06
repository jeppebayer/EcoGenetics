#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 6:00:00

# Sort alignment with regard to QNAME and write output
samtools sort -@ 7 -n -O BAM \
-T people/Jeppe_Bayer/steps/temp/ \
-o people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.namesort.bam \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.raw.bam

exit 0