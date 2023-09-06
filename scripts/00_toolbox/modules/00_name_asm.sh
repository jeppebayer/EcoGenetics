#!/bin/bash

# ----------------- Configuration ----------------------------------------

masking=true
annotation=true

# ------------------ Flag Processing --------------------------------------

while getopts 'f:ma' OPTION; do
    case "$OPTION" in
        f)
            if [ -f "$OPTARG" ]; then
                targetfile="$OPTARG"
            else
                echo "$OPTARG is not a file"
            fi
            ;;
        m)
            masking=false
            ;;
        a)
            annotation=false
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done

# ------------------ Main -------------------------------------------------

[ "$masking" = false ] && mask="_nomask"

[ "$annotation" = false ] && ann="_noann"

if [[ "$targetfile" == */04_assembly/* ]]; then
    path="$targetfile"
    count=0
    dir="$(basename "$targetfile")"
    while [ ! "$dir" == "04_assembly" ]; do
        path="$(dirname "$path")"
        ((count++))
        dir="$(basename "$path")"
    done
    species_name="$targetfile"
    for (( i=1; i<=count-1; i++)); do
        species_name="$(dirname "$species_name")"
    done
    species_name="$(basename "$species_name")"
    genus=${species_name%_*}; genus=${genus::3}; genus=${genus^}
    species=${species_name#*_}; species=${species::3}; species=${species^}
    speciesabbr="$genus""$species"   
fi

referencedir=/faststorage/project/EcoGenetics/BACKUP/reference_genomes/"$species_name"
# Creates specified directory if it doesn't exist
[ -d "$referencedir" ] || mkdir -m 775 "$referencedir"

cp "$targetfile" "$referencedir"/EG_"$speciesabbr"_"$(date +"%d%m%Y")"_genomic"$mask""$ann".fna

exit 0