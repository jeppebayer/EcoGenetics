#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
dataprep=$6 # Path to script location
algo=$7 # Chosen algorithm
age=$8 # Identified category age of sample

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

# If sample is identified as historic uses custom value for -l (seedLen), -n (maxDiff), and -o (maxGapO)
if [ "$age" == " historic" ]; then
    # Gets Suffix Array coordinates for input reads
    bwa "$algo" -t "$cpus" \
    -l 16500 \
    -n 0.01 \
    -o 2 \
    -f "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
    "${RG%.*}" \
    "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated

    bwa "$algo" -t "$cpus" \
    -l 16500 \
    -n 0.01 \
    -o 2 \
    -f "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai \
    "${RG%.*}" \
    "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated

    # Generate alignments
    bwa sampe \
    -r "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
    "${RG%.*}" \
    "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
    "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai \
    "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
    "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
    | \

    # Sort with regards to QNAME and convert to bam format
    samtools sort -@ "$((cpus - 1))" -n -O BAM \
    -T "$WD"/temp/ \
    -o "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
    -

    # File removal
    # rm -f \
    # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
    # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
    # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
    # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai

    exit 0
else
    if [ "$algo" == "mem" ]; then
        # Align sample to reference genome
        bwa "$algo" -t "$cpus" \
        -R "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
        "${RG%.*}" \
        "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
        "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
        | \

        # Sort with regards to QNAME and convert to bam format
        samtools sort -@ "$((cpus - 1))" -n -O BAM \
        -T "$WD"/temp/ \
        -o "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
        -

        # File removal
        # rm -f \
        # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
        # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated

        exit 0
    elif [ "$algo" == "aln" ]; then
        # Gets Suffix Array coordinates for input reads
        bwa "$algo" -t "$cpus" \
        -l 16500 \
        -n 0.01 \
        -o 2 \
        -f "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
        "${RG%.*}" \
        "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated

        bwa "$algo" -t "$cpus" \
        -l 16500 \
        -n 0.01 \
        -o 2 \
        -f "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai \
        "${RG%.*}" \
        "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated

        # Generate alignments
        bwa sampe \
        -r "@RG\tID:$(basename "$sample")\tSM:$(basename "$sample")" \
        "${RG%.*}" \
        "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
        "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai \
        "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
        "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
        | \

        # Sort with regards to QNAME and convert to bam format
        samtools sort -@ "$((cpus - 1))" -n -O BAM \
        -T "$WD"/temp/ \
        -o "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
        -

        # File removal
        # rm -f \
        # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
        # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
        # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.sai \
        # "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.sai

        exit 0
    fi
fi