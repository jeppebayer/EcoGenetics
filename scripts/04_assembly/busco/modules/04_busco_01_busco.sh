#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

targetfile="$1"
WD="$2"
currentuser="$3"
dataset="$4"

if [ "$USER" == "$currentuser" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate busco
fi

name="$(basename "$targetfile")"
name=${name%.*}

busco \
-f \
-i "$targetfile" \
-m genome \
-o busco_"$name" \
--out_path "$WD" \
--download_path /faststorage/project/EcoGenetics/BACKUP/database/busco/busco_downloads \
-l /faststorage/project/EcoGenetics/BACKUP/database/busco/busco_downloads/lineages/"$dataset"

exit 0