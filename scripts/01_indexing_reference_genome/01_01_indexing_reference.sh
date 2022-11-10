#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 2:00:00

# Species specific reference genome, abosolute path (reference genome in FASTA format)
RG="BACKUP/reference_genomes/Pogonognathellus_flavescens/GCA_019776165.1_ASM1977616v1_genomic.fna"

# Indexing reference genome for aligning
bwa index \
-p "${RG%.*}" \
"$RG" ; \

# Adding fai index to reference
samtools faidx \
-o "${RG%.*}".fai \
"$RG"

exit 0