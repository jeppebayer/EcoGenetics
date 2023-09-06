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
    source activate genome_assembly

fi

split_fa \
"$asmfile" \
> "$WD"/"$speciesabbr".asm.split.fasta

exit 0