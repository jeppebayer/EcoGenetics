#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 30:00:00

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ 7 \
-F 0x404 -f 0x2 -q 20 \
-o people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.filtered.bam \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.markdup.bam && \

# Creates files listing number of reads in markdup file on first line
# number of reads in filtered file on second line
# and % of remaining reads on the third line
samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/QC_filter/markdup_to_filtered.readchange && \
samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.filtered.bam \
>> people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/QC_filter/markdup_to_filtered.readchange && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/QC_filter/markdup_to_filtered.readchange \
>> people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/QC_filter/markdup_to_filtered.readchange

exit 0