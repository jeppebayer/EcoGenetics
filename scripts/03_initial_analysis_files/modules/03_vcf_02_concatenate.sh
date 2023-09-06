#!/bin/bash

WD="$1" # Working directory
temp="$2" # Temp directory
sample_name="$3" # Name of current sample
bestn="$4" # Best n alleles value

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

nname=""
if [ $((bestn)) -gt 0 ]; then
    nname="_n=$bestn"
fi

filelist="$temp"/"$sample_name"_filelist.txt
echo -n "" > "$filelist"

for vcf in "$temp"/"$sample_name"_*"$nname"_complete.vcf; do
    echo "$vcf" >> "$filelist"
done

bcftools concat \
--threads 8 \
-f "$filelist" \
-o "$WD"/"$sample_name""$nname"_complete.vcf \
-O v

exit 0