#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

WD="$1"
currentuser="$2"

if [ "$USER" == "$currentuser" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

calcuts \
"$WD"/PB.stat \
> "$WD"/cutoffs

exit 0