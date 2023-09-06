#!/bin/bash

reference_genome="$1"
repmask="$2"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate repeatmasking
fi

species_name="$(basename "$(dirname "$reference_genome")")"

# Creates semi-standard species specific ID code. Ex. turns Aglais_urticae into aglUrt
genus="${species_name%_*}"
genus="${genus::2}"
species="${species_name#*_}"
species="${species::2}"
species_code="$genus""$species"

bedtools \
-soft \
-fi "$reference_genome" \
-bed "$repmask"/"$(basename "$reference_genome")".full_mask.gff \
-fo "$repmask"/"$species_code".full_mask.soft.fna

exit 0