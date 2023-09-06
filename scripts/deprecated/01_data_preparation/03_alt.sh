#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 8:00:00

# Creates bai index for alignment
samtools index -@ 7 -b \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.markdup.bam.bai ; \

# Creates flagstat file for alignment 
samtools flagstat -@ 7 \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/QC_sort/E100049659_L01_ANIcnqkR228559-607.markdup.flagstat ; \

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.markdup.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/QC_sort/E100049659_L01_ANIcnqkR228559-607.markdup.idxstats ; \

# Creates coverage file for alignment
samtools coverage \
-o people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/QC_sort/E100049659_L01_ANIcnqkR228559-607.markdup.coverage \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.markdup.bam

exit 0