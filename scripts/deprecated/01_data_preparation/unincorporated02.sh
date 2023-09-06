#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 8
#SBATCH --time 4:00:00

samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/comparison_of_reads.filter && \
samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam \
>> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/comparison_of_reads.filter && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/comparison_of_reads.filter \
>> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/comparison_of_reads.filter