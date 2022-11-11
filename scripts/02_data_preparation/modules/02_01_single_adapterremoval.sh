#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH --cpus-per-task 8
#SBATCH --time 12:00:00

# Reference genome
RG=$1

# Species directory
SD=$2

# Working directory
WD=$3

# Sample directory
sample=$4

for R in "$sample"/*.fq.gz ; do
    AdapterRemoval \
    --threads 8 \
    --file1 "$R" \
    --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
    --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
    --minquality 25 \
    --minlength 20 \
    --basename "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed \
    --trimns \
    --trimqualities
done ;

# File removal
rm -f \
"$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.discarded \
"$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.settings

exit 0