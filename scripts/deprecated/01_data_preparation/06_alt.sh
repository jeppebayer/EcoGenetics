#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 8:00:00

# Creates flagstat file for alignment 
samtools flagstat -@ 7 \
people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam \
> people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.flagstat ; \

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam \
> people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.idxstats ; \

# Creates coverage file for alignment
samtools coverage \
-o people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_cincta/NYS-F/QC_filter/E100050923_L01_ANIcnqkR228563-642.filtered.coverage \
people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam

exit 0