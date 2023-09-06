#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Creates bai index for filtered alignment
# samtools index -@ 7 -b \
# people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.filtered.bam \
# > people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.filtered.bam.bai && \

# Make mpileup file
samtools mpileup \
-C 50 \
-o people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.mpileup \
-f BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna \
people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.filtered.bam

exit 0