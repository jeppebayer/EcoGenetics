#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 8
#SBATCH --time 4:00:00

samtools flagstat -@ 7 \
people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.flagstat \
samtools idxstats -@ 7 \
people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.idxstats \
samtools coverage \
people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.coverage