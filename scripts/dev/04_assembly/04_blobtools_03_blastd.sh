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

diamond blastx \
--query "$asm_file" \
--max-target-seqs 1 \
--sensitive \
--threads 16 \
--db /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/uniprot/uniprot_trembl.dmnd \
--evalue 1e-25 \
--outfmt 6 \
--out "$WD"/"$speciesabbr"_diamondblast.out

# taxify results
blobtools taxify \
-f "$WD"/"$speciesabbr"_diamondblast.out \
-m uniprot_ref_proteomes.taxids \
-s 0 \
-t 2 \
-o "$WD"/

exit 0