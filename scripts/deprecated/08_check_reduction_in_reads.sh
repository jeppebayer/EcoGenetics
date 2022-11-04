#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 8
#SBATCH --time 4:00:00

samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/comparison_of_reads.filtered && \
samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/E100050923_L01_ANIcnqkR228563-642.bam \
>> people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/comparison_of_reads.filtered && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/comparison_of_reads.filtered \
>> people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/comparison_of_reads.filtered