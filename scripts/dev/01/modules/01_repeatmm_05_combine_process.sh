#!/bin/bash

reference_genome="$1"
repbaserun="$2"
repmodrun="$3"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate repeatmasking
fi

cat \
"$repbaserun"/"$(basename "$reference_genome")".cat.gz \
"$repmodrun"/"$(basename "$reference_genome")".cat.gz \
> "$(basename "$reference_genome")".full_mask.cat.gz

cat \
"$repbaserun"/"$(basename "$reference_genome")".out \
"$repmodrun"/"$(basename "$reference_genome")".out \
> "$(basename "$reference_genome")".full_mask.out

ProcessRepeats \
-gff \
-dir /faststorage/project/EcoGenetics/BACKUP/database/giri_repbase/RepBAseCustom15.02.fast/Arthropoda.rep \
"$(basename "$reference_genome")".full_mask.cat.gz

exit 0