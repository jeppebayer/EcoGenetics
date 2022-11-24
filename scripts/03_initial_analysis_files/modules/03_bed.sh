#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
dataprep=$6 # Path to script location

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

# Creates .bed file of intergenic regions
bedtools complement \
-i "$WD"/temp/"$(basename "${RG%.*}")"_gene_sorted.bed \
-g "$WD"/temp/"$(basename "${RG%.*}")"_genome_sorted.txt \
> "$WD"/temp/"$(basename "${RG%.*}")"_intergenic.bed # Check if the intergenic bed file should be saved in ref genome folder

# File removal
rm -f \
"$WD"/temp/"$(basename "${RG%.*}")"_gene_sorted.bed \
"$WD"/temp/"$(basename "${RG%.*}")"_genome_sorted.txt

exit 0