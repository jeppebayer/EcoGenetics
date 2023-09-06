#!/bin/bash

WD="$1" # Working directory
sample="$2" # Sample directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Merge R1R2 and collapsed file to one name-sorted bam file
samtools merge -@ 7 \
-o "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
-c -p -n \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_collapsed_aligned.bam

# File removal
rm -f \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_collapsed_aligned.bam

exit 0