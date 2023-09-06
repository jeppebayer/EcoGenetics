#!/bin/bash

reference_genome="$1"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate repeatmasking
fi

species_name="$(basename "$(dirname "$reference_genome")")"

# Creates semi-standard species specific ID code. Ex. turns Aglais_urticae into AglUrt
genus="${species_name%_*}"
genus="${genus::3}"
genus="${genus^}"
species="${species_name#*_}"
species="${species::3}"
species="${species^}"
species_code="$genus""$species"1

# Actually runs RepeatModeler
RepeatModeler \
-database ./RM_DB_"$species_name"/"$species_name" \
-LTRStruct \
-pa 8
# -threads 32 \

# Adds species code to outputted fasta file
cat ./*/"$species_name"-families.fa | \
seqkit fx2tab | \
awk -v species_code="$species_code" '{ print species_code"_"$0 }' | \
seqkit tab2fx \
> "$species_name"-families.prefix.fa

# Split fasta file into classified and unclassified elements
cat "$species_name"-families.prefix.fa | \
seqkit fx2tab | \
rg -j 10 -v "Unknown" | \
seqkit tab2fx \
> "$species_name"-families.prefix.fa.known

cat "$species_name"-families.prefix.fa | \
seqkit fx2tab | \
rg -j 10 "Unknown" | \
seqkit tab2fx \
> "$species_name"-families.prefix.fa.unknown

exit 0