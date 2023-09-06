#!/bin/bash

WD="$1" # Working directory
sample="$2" # Sample directory
temp="$3" # Temporary file directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Concatenation of R1 files
R1_1=
for R1_2 in "$sample"/*_1.fq.gz ; do
    if [ "$R1_1" ] && [ ! -f "$temp"/"$(basename "$sample")"_R1.fq.gz ]; then
        cat "$R1_1" "$R1_2" > "$temp"/"$(basename "$sample")"_R1.fq.gz
    elif [ -f "$temp"/"$(basename "$sample")"_R1.fq.gz ]; then
        cat "$R1_2" >> "$temp"/"$(basename "$sample")"_R1.fq.gz
    else
        R1_1=$R1_2
    fi
done

AdapterRemoval \
--threads 8 \
--file1 "$temp"/"$(basename "$sample")"_R1.fq.gz \
--adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
--adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
--minquality 25 \
--minlength 20 \
--basename "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed \
--trimns \
--trimqualities

# File removal
rm -f \
"$temp"/"$(basename "$sample")"_R1.fq.gz \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.discarded \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.settings

exit 0