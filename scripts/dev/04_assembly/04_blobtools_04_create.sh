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

blobtools create \
--nodes /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/ncbi_taxonomy/nodes.dmp \
--names /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/ncbi_taxonomy/names.dmp \
-i "$asm_file" \
-t "$WD"/"$speciesabbr"_diamondblast.taxified.out \
-t "$WD"/"$speciesabbr"_ncbimegablast.out \
-x bestsumorder \
-c "$WD"/"$speciesabbr".bam.cov \
-o "$WD"/"$speciesabbr"_blobtools_db

exit 0