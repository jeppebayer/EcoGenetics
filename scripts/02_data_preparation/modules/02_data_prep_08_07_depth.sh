#!/bin/bash

SD="$1" # Species directory
WD="$2" # Working directory
sample="$3" # Sample directory
script_path="$4" # Script path

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

samtools depth \
-@ 7 \
-a \
-o "$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.depth \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam

# Creates histgram of read depth
python "$script_path"/modules/02_data_prep_08_07_plotdepth.py \
"$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.depth \
"$sample"

rm "$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.depth

exit 0