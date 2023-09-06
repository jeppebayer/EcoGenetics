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

blobtools create \
--nodes /faststorage/project/EcoGenetics/BACKUP/database/ncbi_taxonomy/nodes.dmp \
--names /faststorage/project/EcoGenetics/BACKUP/database/ncbi_taxonomy/names.dmp \
-i "$asmfile" \
-t "$WD"/"$speciesabbr"_diamondblast.taxified.out \
-t "$WD"/"$speciesabbr"_ncbimegablast.out \
-x bestsumorder \
-c "$WD"/"$speciesabbr".bam.cov \
-o "$WD"/"$speciesabbr"_blobtools_db

blobtools plot \
-i "$WD"/"$speciesabbr"_blobtools_db.blobDB.json \
-o "$WD"/ \
-x bestsumorder

blobtools view \
-i "$WD"/"$speciesabbr"_blobtools_db.blobDB.json \
-o "$WD"/ \
-x bestsumorder \
-r all

exit 0