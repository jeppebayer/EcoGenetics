#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

asmfile="$1"
WD="$2"
speciesabbr="$3"
currentuser="$4"

if [ "$USER" == "$currentuser" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools

fi

blobtools map2cov \
-i "$asmfile" \
-b "$WD"/"$speciesabbr".bam \
-o "$WD"/

exit 0