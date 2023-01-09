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

blastn \
-task megablast \
-query "$asm_file" \
-db /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/blastdb/nt \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 10 \
-max_hsps 1 \
-num_threads 16 \
-evalue 1e-25 \
-out "$WD"/"$speciesabbr"_ncbimegablast.out

exit 0