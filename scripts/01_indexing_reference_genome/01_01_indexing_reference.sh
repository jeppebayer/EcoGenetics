#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 2:00:00

# Indexing reference genome for aligning
bwa index \
-p BACKUP/reference_genomes/Pogonognathellus_flavescens/GCA_019776165.1_ASM1977616v1_genomic \
BACKUP/reference_genomes/Pogonognathellus_flavescens/GCA_019776165.1_ASM1977616v1_genomic.fna ; \

# Adding fai index to reference
samtools faidx \
-o BACKUP/reference_genomes/Pogonognathellus_flavescens/GCA_019776165.1_ASM1977616v1_genomic.fna.fai \
BACKUP/reference_genomes/Pogonognathellus_flavescens/GCA_019776165.1_ASM1977616v1_genomic.fna

exit 0