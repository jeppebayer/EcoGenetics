#!/bin/bash

genome="$1"
WD="$(dirname "$1")"

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate genome_assembly
fi


# Creates files with all sequence lengths sorted longest to shortest
bioawk \
-v OFS='\t' \
-c fastx \
'{ print $name, length($seq) }' \
"$genome" \
| sort -k1,1 -k2,2nr \
> "$WD"/sequence_lengths

# Calculates total assembly length
total=$(awk '{Total=Total+$2} END{print Total}' "$WD"/sequence_lengths)
awk '{Total=Total+$2} END{print "Total assembly length is: " Total}' "$WD"/sequence_lengths > "$WD"/asm_lenght_"$total"