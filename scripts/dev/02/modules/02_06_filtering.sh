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

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ "$((cpus - 1))" \
-F 3844 -q 20 \
-o "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ "$((cpus - 1))" -b \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam.bai

exit 0