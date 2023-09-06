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

for R in "$sample"/*.fq.gz ; do
    AdapterRemoval \
    --threads 8 \
    --file1 "$R" \
    --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
    --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
    --minquality 25 \
    --minlength 20 \
    --basename "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed \
    --trimns \
    --trimqualities
done ;

# File removal
# rm -f \
# "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.discarded \
# "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.settings

exit 0