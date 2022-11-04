#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Creates bai index for alignment
samtools index -@ 7 -b \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/E100050923_L01_ANIcnqkR228563-642.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/E100050923_L01_ANIcnqkR228563-642.markdup.bam.bai ; \

# Creates flagstat file for alignment 
samtools flagstat -@ 7 \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/E100050923_L01_ANIcnqkR228563-642.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.flagstat ; \

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/E100050923_L01_ANIcnqkR228563-642.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.idxstats ; \

# Creates coverage file for alignment
samtools coverage \
-o people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.coverage \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/E100050923_L01_ANIcnqkR228563-642.markdup.bam

exit 0