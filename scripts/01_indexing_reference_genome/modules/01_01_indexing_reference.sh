#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

# Reference genome
RG=$1

# Indexing reference genome for aligning
bwa index \
-p "${RG%.*}" \
"$RG" ; \

# Adding fai index to reference
samtools faidx \
-o "$RG".fai \
"$RG"

exit 0