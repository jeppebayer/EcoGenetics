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

# Creates text file listing number of reads in markdup file on first line
# number of reads in filtered file on second line and % remaining reads on the third line
samtools view -@ 7 -c \
"$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange && \
samtools view -@ 7 -c \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
>> "$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
"$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange \
>> "$WD"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange

# File removal
# rm -f \
# "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

exit 0