#!/bin/bash

WD="$1" # Working directory
sample="$2" # Sample directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

R1=

for R2 in "$sample"/*.fq.gz ; do
    if [ "$R1" ] ; then
        # 2nd entry
        AdapterRemoval \
        --threads 8 \
        --file1 "$R1" \
        --file2 "$R2" \
        --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        --minquality 25 \
        --minlength 20 \
        --basename "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed \
        --trimns \
        --trimqualities \
        --collapse

    else
        # First entry
        R1=$R2
    fi
done

# File removal
rm -f \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.discarded \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.settings \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.singleton.truncated

exit 0