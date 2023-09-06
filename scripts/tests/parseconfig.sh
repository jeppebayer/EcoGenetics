#!/bin/bash

config_file="$(readlink -f "$1")"

# Use '\t' as field separater
# Skip lines starting with '#', a blank space or '>'
# Substitute all spaces with '_'
settings=($(awk \
    -F "\t" \
    '/^[^#^\ ^>]/ {
        print gensub(/\ /, "_", "g")
    }' \
    "$config_file"))

# order=($(awk \
#     -F "\t" \
#     '/^[>]/ {
#         print tolower(gensub(/\ /, "_", "g", gensub(/[>:]/, "", "g")))
#     }' \
#     "$config_file"))

species_name="${settings[0]}"
genus="${species_name%_*}"; genus="${genus::3}"; genus="${genus^}"
if [ "${genus: -1}" == "." ]; then
    genus=${genus::-1}
fi
species="${species_name#*_}"; species="${species::3}"; species="${species^}"
if [ "${species: -1}" == "." ]; then
    species=${species::-1}
fi
species_code="$genus""$species"

echo "$species_code"

reference_genome="$(readlink -f "${settings[1]}")"

echo "$reference_genome"

reseq_dir="$(readlink -f "${settings[2]}")"

echo "$reseq_dir"

working_dir="$(readlink -f "${settings[3]}")"

echo "$working_dir"



exit 0