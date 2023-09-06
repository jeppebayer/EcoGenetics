#!/bin/bash

SD="$1" # Species directory
WD="$2" # Working directory
sample="$3" # Sample directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ 7 \
-F 3844 -q 20 \
-o "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ 7 -b \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam.bai

exit 0