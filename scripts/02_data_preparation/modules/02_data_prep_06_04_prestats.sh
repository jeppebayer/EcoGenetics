#!/bin/bash

WD="$1" # Working directory
sample="$2" # Sample directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Creates stat file for alignment
samtools stats -@ 7 \
-c 1,1000,1 \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.stats

exit 0