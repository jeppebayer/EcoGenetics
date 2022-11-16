#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Reference genome
RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Maniola_jurtina/GCF_905333055.1_ilManJurt1.1_genomic.fna"

# Species directory
SD="/home/jepe/EcoGenetics/BACKUP/museomics/Maniola_jurtina"

# Working directory
WD=$3

# Sample directory
sample="/home/jepe/EcoGenetics/BACKUP/museomics/Maniola_jurtina/MaJu_01_J_2022"

export DISPLAY=:0

# First checks whether a .gff file for the reference genome is available
for file in "$(dirname "$RG")"/*.gff; do
    if [ -e "$file" ]; then
        
        # .gff file is available
        gff=$file

        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -gff "$gff" \
        -outdir /home/jepe/EcoGenetics/people/Jeppe_Bayer/data/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=16G
        exit 0

    else

        # .gff file is not available
        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -outdir /home/jepe/EcoGenetics/people/Jeppe_Bayer/data/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=16G
        exit 0
    fi
done
