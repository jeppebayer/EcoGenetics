#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools

fi

asm_file=$1
hifi_reads=$2
WD=$3
temp=$4
speciesabbr=$5

blobtools plot \
-i "$WD"/"$speciesabbr"_blobtools_db.blobDB.json \
-o "$WD"/"$speciesabbr"

exit 0