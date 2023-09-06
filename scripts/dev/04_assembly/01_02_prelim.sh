#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --time=60
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Entomobrya_nicoleti/blobtools/out/prelim-%j.out \

asmfile="$1"
WD="$2"
speciesabbr="$3"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools
fi

# blobtools create \
# --nodes /faststorage/project/EcoGenetics/BACKUP/database/ncbi_taxonomy/nodes.dmp \
# --names /faststorage/project/EcoGenetics/BACKUP/database/ncbi_taxonomy/names.dmp \
# -i "$asmfile" \
# -t "$WD"/"$speciesabbr"_ncbimegablast.out \
# -c "$WD"/"$speciesabbr".bam.cov \
# -o "$WD"/"$speciesabbr"_blobtools_db

# blobtools plot \
# -i "$WD"/"$speciesabbr"_blobtools_db.blobDB.json \
# -o "$WD"/

# blobtools view \
# -i "$WD"/"$speciesabbr"_blobtools_db.blobDB.json \
# -o "$WD"/

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