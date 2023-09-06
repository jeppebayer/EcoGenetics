#!/bin/bash

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
dataprep=$6 # Path to script location
algo=$7 # Chosen algorithm

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ 7 \
-f 4 \
-o "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_unmapped.bam \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ 7 -b \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_unmapped.bam \
> "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_unmapped.bam.bai

exit 0