#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 8
#SBATCH --time 10:00:00

samtools fixmate -@ 7 -m -O BAM \
people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.fixmate.bam && \
samtools markdup -@ 7 -s \
-f people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/stats.fixmate \
-T people/Jeppe_Bayer/steps/temp/temp. \
people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.fixmate.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/E100050923_L01_ANIcnqkR228563-642.markdup.bam