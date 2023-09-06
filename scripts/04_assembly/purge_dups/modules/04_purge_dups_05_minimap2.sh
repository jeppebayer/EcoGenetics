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

minimap2 \
-x asm5 \
-t 15 \
-DP \
"$WD"/"$speciesabbr".asm.split.fasta \
"$WD"/"$speciesabbr".asm.split.fasta \
| gzip \
-c \
- \
> "$WD"/"$speciesabbr"_alignment.self.paf.gz

exit 0