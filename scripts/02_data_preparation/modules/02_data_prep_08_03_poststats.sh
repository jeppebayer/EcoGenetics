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

# Creates coverage file for alignment
samtools coverage \
-o "$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.coverage \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam

exit 0