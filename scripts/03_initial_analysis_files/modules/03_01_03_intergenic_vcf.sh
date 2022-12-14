#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH --cpus-per-task 15
#SBATCH --time 30:00:00
#SBATCH --output=/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps/03_init_analysis_files/intergenic_vcf-%j.out

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
data="03_init_analysis_files"
SD="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta"

# Creates .vcf files for sample
freebayes-parallel \
<(fasta_generate_regions.py "$RG" --chunks 100) 15 \
-t "$(dirname "$RG")"/"$(basename "${RG%.*}")"_intergenic.bed \
-p 100 \
-f "$RG" \
-v "$WD"/"$data"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_intergenic.vcf \
--pooled-discrete \
"$sample"/"$(basename "$sample")"_filtered.bam

# Creates .vcf files for sample
# freebayes-parallel \
# <(fasta_generate_regions.py "$RG" --chunks "$cpus") "$cpus" \
# -t "$(dirname "$RG")"/"$(basename "${RG%.*}")"_intergenic.bed \
# -p 100 \
# -f "$RG" \
# -v "$WD"/"$data"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_intergenic.vcf \
# --pooled-discrete \
# "$sample"/"$(basename "$sample")"_filtered.bam

# Annotates .vcf file
# snpEff ann > "$WD"/03_initial_analysis_files/"$(basename "$sample")"_intergenic.vcf

# freebayes -f "$RG" \
# -p 100 \
# -f "$RG" \
# -v "$sample"/"$(basename "$sample")".vcf \
# --pooled-discrete \
# "$sample"/"$(basename "$sample")"_filtered.bam

exit 0