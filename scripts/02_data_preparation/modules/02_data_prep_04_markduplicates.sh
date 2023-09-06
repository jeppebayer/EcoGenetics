#!/bin/bash

WD="$1" # Working directory
sample="$2" # Sample directory
temp="$3" # Temp directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Add fixmate tag to alignment. Can only be done on name sorted alignment
samtools fixmate -@ 7 -m -O BAM \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
- | \

# Position sort alignment and pipe output
samtools sort -@ 7 -O BAM \
-T "$temp"/ \
- | \

# Mark duplicates and save output as .coordinatesort.bam. Also outputs some realted statistics
samtools markdup -@ 7 -s \
-f "$WD"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.markdupstats \
-T "$temp"/ \
- \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

# File removal
rm -f \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam

exit 0