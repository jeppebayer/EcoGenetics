#!/bin/bash

WD="$1" # Working directory
sample="$2" # Sample directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

for R in "$sample"/*.fq.gz ; do
    AdapterRemoval \
    --threads 8 \
    --file1 "$R" \
    --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
    --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
    --minquality 25 \
    --minlength 20 \
    --basename "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed \
    --trimns \
    --trimqualities
done

# File removal
rm -f \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.discarded \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.settings

exit 0