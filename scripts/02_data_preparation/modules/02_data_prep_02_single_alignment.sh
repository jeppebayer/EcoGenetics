#!/bin/bash

RG="$1" # Reference genome
WD="$2" # Working directory
sample="$3" # Sample directory
temp="$4" # Temp directory
algo="$5" # Chosen algorithm
age="$6" # Identified category age of sample

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Checks naming format of reference index and tires to match
if [ -e "$RG".ann ]; then
    idx="$RG"
else
    idx="${RG%.*}"
fi

# If sample is identified as historic uses custom value for -l (seedLen), -n (maxDiff), and -o (maxGapO)
if [ "$age" == " historic" ]; then
    # Gets Suffix Array coordinates for input reads
    bwa "$algo" -t 8 \
    -l 16500 \
    -n 0.01 \
    -o 2 \
    "$idx" \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.truncated \
    | \

    # Generate alignments
    bwa samse \
    -r "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
    "$idx" \
    - \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.truncated \
    | \

    # Sort with regards to QNAME and convert to bam format
    samtools sort -@ 7 -n -O BAM \
    -T "$temp"/ \
    -o "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
    -

    # File removal
    rm -f \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.truncated

    exit 0
else
    if [ "$algo" == "mem" ]; then
        # Align sample to reference genome
        bwa "$algo" -t 8 \
        -R "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
        "$idx" \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.truncated \
        | \

        # Sort with regards to QNAME and convert to bam format
        samtools sort -@ 7 -n -O BAM \
        -T "$WD"/temp/ \
        -o "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
        -

        # File removal
        rm -f \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.truncated

        exit 0
    elif [ "$algo" == "aln" ]; then
        # Gets Suffix Array coordinates for input reads
        bwa "$algo" -t 8 \
        "$idx" \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.truncated \
        | \

        # Generate alignments
        bwa samse \
        -r "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
        "$idx" \
        - \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.truncated \
        | \

        # Sort with regards to QNAME and convert to bam format
        samtools sort -@ 7 -n -O BAM \
        -T "$WD"/temp/ \
        -o "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
        -

        # File removal
        rm -f \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.truncated

        exit 0
    fi
fi
