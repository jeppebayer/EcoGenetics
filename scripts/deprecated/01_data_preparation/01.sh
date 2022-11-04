#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Indexing reference genome for aligning and for later conversion of data to vcf format

bwa index \
-p /home/jepe/EcoGenetics/BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic \
-a bwtsw \
/home/jepe/EcoGenetics/BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna ; \
samtools faidx \
-o BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna.fai \
BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna

exit 0