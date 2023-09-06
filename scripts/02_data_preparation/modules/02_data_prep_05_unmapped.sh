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

# Removes all mapped reads
samtools view -b -@ 7 \
-f 4 \
-o "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_unmapped.bam \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ 7 -b \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_unmapped.bam \
> "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_unmapped.bam.bai

exit 0