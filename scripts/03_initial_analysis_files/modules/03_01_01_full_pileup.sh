#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 30:00:00
#SBATCH --output=/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps/03_initial_analysis_files/03_01_mpileup-%j.out

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
data=$6 # Path to data location in WD

RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna"
sample="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F"

# Make pileup file
samtools mpileup \
-C 0 \
-o "$sample"/"$(basename "$sample")".pileup \
-f "$RG" \
"$sample"/"$(basename "$sample")"_filtered.bam

exit 0