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

# Merge R1R2 and collapsed file to one name-sorted bam file
samtools merge -@ "$((cpus - 1))" \
-o "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
-c -p -n \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_collapsed_aligned.bam

# File removal
# rm -f \
# "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
# "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_collapsed_aligned.bam

exit 0