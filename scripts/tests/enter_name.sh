#!/bin/bash

read -rp "Enter species name: " species_name

if [[ "$species_name" == *\ * ]]; then
    genus="${species_name%\ *}"; genus="${genus::3}"; genus="${genus^}"
    if [ "${genus: -1}" == "." ]; then
        genus=${genus::-1}
    fi
    species="${species_name#*\ }"; species="${species::3}"; species="${species^}"
    if [ "${species: -1}" == "." ]; then
        species=${species::-1}
    fi
    species_code="$genus""$species"
    echo "$species_code"
elif [[ "$species_name" == *_* ]]; then
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
else
    echo "Species name has to be provided as '<genus> <species>' | '<genus>_<species>"
    exit 1
fi

# echo "$(dirname "$0")"

exit 0