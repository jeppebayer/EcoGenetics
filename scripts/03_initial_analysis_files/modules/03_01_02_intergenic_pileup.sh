#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 30:00:00
#SBATCH --dependency=afterany:10558521
#SBATCH --output=/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps/03_initial_analysis_files/03_intergenic_pileup-%j.out

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
data=$6 # Path to data location in WD

RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna"
sample="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F"
WD="/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps"

if [ ! -e "$(dirname "$RG")"/"$(basename "${RG%.*}")"_intergenic.bed ]; then

    # Creates sorted .bed files containing all start and end position of all gene regions from .gff file
    grep "gene" "$(dirname "$RG")"/"$(basename "${RG%.*}")".gff \
    | awk '{OFS="\t"}{print $1, $4, $5}' \
    | sort -k1,1 -k2,2n \
    > "$WD"/temp/"$(basename "${RG%.*}")"_gene_sorted.bed

    # Creates sorted file with the length of all scaffolds
    bioawk \
    -v OFS='\t' \
    -c fastx \
    '{ print $name, length($seq) }' \
    "$RG" \
    | sort -k1,1 -k2,2n \
    > "$WD"/temp/"$(basename "${RG%.*}")"_genome_sorted.txt

    # Creates .bed file of intergenic regions and places in in species ref genome directory
    bedtools complement \
    -i "$WD"/temp/"$(basename "${RG%.*}")"_gene_sorted.bed \
    -g "$WD"/temp/"$(basename "${RG%.*}")"_genome_sorted.txt \
    > "$(dirname "$RG")"/"$(basename "${RG%.*}")"_intergenic.bed

    # File removal
    rm -f \
    "$WD"/temp/"$(basename "${RG%.*}")"_gene_sorted.bed \
    "$WD"/temp/"$(basename "${RG%.*}")"_genome_sorted.txt

fi

# Make intergenic pileup file
samtools mpileup \
-l "$(dirname "$RG")"/"$(basename "${RG%.*}")"_intergenic.bed \
-C 0 \
-o "$WD"/"$data"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_intergenic.pileup \
-f "$RG" \
"$sample"/"$(basename "$sample")"_filtered.bam

exit 0