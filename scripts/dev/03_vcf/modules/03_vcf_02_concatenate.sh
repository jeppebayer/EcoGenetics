#!/bin/bash

WD="$1" # Working directory
temp="$2" # Temp directory
sample_name="$3" # Name of current sample

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

filelist="$temp"/"$sample_name"_filelist.txt
touch "$filelist"

for vcf in "$temp"/"$sample_name"_*_complete.vcf; do
    echo "$vcf" >> "$filelist"
done

bcftools concat \
--threads 6 \
-f "$filelist" \
-o "$WD"/"$sample_name".vcf \
-O v