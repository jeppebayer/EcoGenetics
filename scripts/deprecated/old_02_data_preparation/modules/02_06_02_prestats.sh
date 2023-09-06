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

# Creates flagstat file for alignment 
samtools flagstat -@ "$((cpus - 1))" \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.flagstat

exit 0