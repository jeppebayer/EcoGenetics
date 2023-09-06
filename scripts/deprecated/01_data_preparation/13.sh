#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 8:00:00

# Creates files listing number of reads in markdup file on first line
# number of reads in filtered file on second line
# and % of remaining reads on the third line
samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/markdup_to_filtered.readchange && \
samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam \
>> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/markdup_to_filtered.readchange && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/markdup_to_filtered.readchange \
>> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/markdup_to_filtered.readchange