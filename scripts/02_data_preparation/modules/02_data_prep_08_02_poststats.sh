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

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.idxstats

exit 0