#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

rawreads="$1"
asmfile="$2"
WD="$3"
speciesabbr="$4"
currentuser="$5"

if [ "$USER" == "$currentuser" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

minimap2 \
-x map-hifi \
-t 15 \
"$asmfile" \
"$rawreads" \
| gzip \
-c \
- \
> "$WD"/"$speciesabbr"_alignment.paf.gz

exit 0