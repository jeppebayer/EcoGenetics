#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Make mpileup file
samtools mpileup \
-C 50 \
-o people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.mpileup \
-f BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna \
people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam

exit 0