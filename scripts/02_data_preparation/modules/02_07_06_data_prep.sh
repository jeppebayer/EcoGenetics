#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Reference genome
# RG=$1

# Species directory
# SD=$2

# Working directory
# WD=$3

# Sample directory
# sample=$4

for file in "$(dirname "$RG")"/*.gff; do
    gff=$file
done

WD="/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps"

SD="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta"

sample="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F"

qualimap bamqc \
-bam "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
-gff "$gff" \
-outdir "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/qualimap/ \
-outfile "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_qualimap.pdf \
-outformat PDF \
--java-mem-size=16G

exit 0