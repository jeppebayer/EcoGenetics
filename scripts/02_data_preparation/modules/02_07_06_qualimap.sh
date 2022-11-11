#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Reference genome
RG=$1

# Species directory
SD=$2

# Working directory
WD=$3

# Sample directory
sample=$4

# First checks whether a .gff file for the reference genome is available
for file in "$(dirname "$RG")"/*.gff; do
    if [ -e "$file" ]; then
        
        # .gff file is available
        gff=$file

        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -gff "$gff" \
        -outdir "$SD"/"$(basename "$sample")"/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=16G
        exit 0

    else

        # .gff file is not available
        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -outdir "$SD"/"$(basename "$sample")"/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=16G
        exit 0
    fi
done
