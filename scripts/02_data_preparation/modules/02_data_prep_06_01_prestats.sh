#!/bin/bash

WD="$1" # Working directory
sample="$2" # Sample directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Creates bai index for alignment
samtools index -@ 7 -b \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam.bai ; \

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.idxstats

# File removal
rm -f \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam.bai

exit 0