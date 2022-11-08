#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 30:00:00

# Make mpileup file
samtools mpileup \
-C 50 \
-o people/Jeppe_Bayer/data/Ocin_NYS-F.pileup \
-f /home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna \
/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps/01_data_preparation/Orchesella_cincta/Ocin_NYS-F/Ocin_NYS-F_filtered.bam
exit 0