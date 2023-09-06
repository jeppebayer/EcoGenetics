#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

WD="$1"
speciesabbr="$2"
currentuser="$3"

if [ "$USER" == "$currentuser" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

pbcstat \
-O "$WD" \
"$WD"/"$speciesabbr"_alignment.paf.gz

exit 0