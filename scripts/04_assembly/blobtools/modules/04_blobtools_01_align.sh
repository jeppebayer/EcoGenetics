#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

rawreads="$1"
asmfile="$2"
WD="$3"
temp="$4"
speciesabbr="$5"
currentuser="$6"

if [ "$USER" == "$currentuser" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools

fi

minimap2 \
-x map-hifi \
-a \
-t 16 \
"$asmfile" \
"$rawreads" \
| samtools sort \
-@ 15 \
-o "$WD"/"$speciesabbr".bam \
-O BAM \
-T "$temp" \
-

samtools index \
-b \
-@ 15 \
"$WD"/"$speciesabbr".bam 

exit 0