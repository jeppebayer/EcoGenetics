#!/bin/bash

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
dataprep=$6 # Path to script location
algo=$7 # Chosen algorithm

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

# Creates text file listing number of reads in markdup file on first line
# number of reads in filtered file on second line and % remaining reads on the third line
samtools view -@ "$((cpus - 1))" -c \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange && \
samtools view -@ "$((cpus - 1))" -c \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
>> "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange \
>> "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange

# File removal
# rm -f \
# "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

exit 0