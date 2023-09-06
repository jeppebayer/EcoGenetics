#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

RG=$1 # Reference genome

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep

fi

# Indexing reference genome for aligning
bwa index \
-p "${RG%.*}" \
"$RG" ; \

# Adding fai index to reference
samtools faidx \
-o "$RG".fai \
"$RG"

exit 0