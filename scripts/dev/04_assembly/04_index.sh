#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

if [ "$USER" == "jepe" ]; then
    source /home/"$USER"/.bashrc
    source activate genome_assembly
fi

# Reference genome
fasta=$1

# Indexing reference genome for aligning
bwa index \
-p "${fasta%.*}" \
"$fasta"

exit 0