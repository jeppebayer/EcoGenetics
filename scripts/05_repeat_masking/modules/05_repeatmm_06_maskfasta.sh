#!/bin/bash

reference_genome="$1"
repmask="$2"
WD="$3"

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
species_code="$genus""$species"

filename="$(basename "$reference_genome")"
filename=${filename%.*}

bedtools maskfasta \
    -soft \
    -fi "$reference_genome" \
    -bed "$repmask"/"$filename".repeats.gff \
    -fo "$WD"/"$species_code".full_mask.soft.fna

# Make BED3 file of repeat regions
awk \
    -F "\t" \
    'BEGIN {OFS = "\t"}
    {if ($0 ~ /^[^#]/)
        {print $1, ($4 - 1), $5}
    }' \
    "$repmask"/"$filename".repeats.gff \
    > "$WD"/"$species_code".repeats.bed

exit 0