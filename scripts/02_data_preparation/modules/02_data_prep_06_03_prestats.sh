#!/bin/bash

WD="$1" # Working directory
sample="$2" # Sample directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Creates coverage file for alignment
samtools coverage \
-o "$WD"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.coverage \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

exit 0