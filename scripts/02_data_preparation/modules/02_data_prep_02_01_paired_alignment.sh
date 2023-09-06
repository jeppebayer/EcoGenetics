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

# If sample is identified as historic uses custom value for -l (seedLen), -n (maxDiff), and -o (maxGapO)
if [ "$age" == " historic" ]; then
    # Gets Suffix Array coordinates for input reads
    bwa "$algo" -t 8 \
    -l 16500 \
    -n 0.01 \
    -o 2 \
    -f "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
    "${RG%.*}" \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated

    bwa "$algo" -t 8 \
    -l 16500 \
    -n 0.01 \
    -o 2 \
    -f "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai \
    "${RG%.*}" \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated

    # Generate alignments
    bwa sampe \
    -r "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
    "${RG%.*}" \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
    | \

    # Sort with regards to QNAME and convert to bam format
    samtools sort -@ 7 -n -O BAM \
    -T "$temp"/ \
    -o "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
    -

    # File removal
    rm -f \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
    "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai

    exit 0
else
    if [ "$algo" == "mem" ]; then
        # Align sample to reference genome
        bwa "$algo" -t 8 \
        -R "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
        "${RG%.*}" \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
        | \

        # Sort with regards to QNAME and convert to bam format
        samtools sort -@ 7 -n -O BAM \
        -T "$temp"/ \
        -o "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
        -

        # File removal
        rm -f \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated

        exit 0
    elif [ "$algo" == "aln" ]; then
        # Gets Suffix Array coordinates for input reads
        bwa "$algo" -t 8 \
        -l 16500 \
        -n 0.01 \
        -o 2 \
        -f "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
        "${RG%.*}" \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated

        bwa "$algo" -t 8 \
        -l 16500 \
        -n 0.01 \
        -o 2 \
        -f "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai \
        "${RG%.*}" \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated

        # Generate alignments
        bwa sampe \
        -r "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
        "${RG%.*}" \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
        | \

        # Sort with regards to QNAME and convert to bam format
        samtools sort -@ 7 -n -O BAM \
        -T "$temp"/ \
        -o "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
        -

        # File removal
        rm -f \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
        "$WD"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai

        exit 0
    fi
fi