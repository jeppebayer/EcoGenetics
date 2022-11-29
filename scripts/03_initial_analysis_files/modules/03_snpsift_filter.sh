#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 30:00:00
#SBATCH --dependency=afterany:10559103
#SBATCH --output=/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps/03_initial_analysis_files/03_snpsift_filter-%j.out

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
data=$6 # Path to data location in WD
n=$7 # Number of parts
script_path=$8 # Path to script location

RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna"
sample="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F"
WD="/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps"

SnpSift filter "ANN[*].EFFECT has 'missense_variant'" \
"$sample"/"$(basename "$sample")"_ann.vcf \
> "$WD"/"$data"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_nonsyn.vcf

exit 0