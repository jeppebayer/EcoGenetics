#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
data=$6 # Path to data location in WD
n=$7 # Number of parts
script_path=$8 # Path to script location

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