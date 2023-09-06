#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --time=8640 \
#SBATCH --mem-per-cpu=15G \
#SBATCH --cpus-per-task=16 \
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Lepidocyrtus_sp/blobtools/out/blastd-%j.out \

asmfile="$1"
WD="$2"
speciesabbr="$3"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools
fi

diamond blastx \
--query "$asmfile" \
--max-target-seqs 1 \
--sensitive \
--threads 16 \
--db /faststorage/project/EcoGenetics/BACKUP/database/uniprot/uniprot_ref_proteomes.dmnd \
--evalue 1e-25 \
--outfmt 6 \
--out "$WD"/"$speciesabbr"_diamondblast.out

# taxify results
blobtools taxify \
-f "$WD"/"$speciesabbr"_diamondblast.out \
-m /faststorage/project/EcoGenetics/BACKUP/database/uniprot/uniprot_ref_proteomes.taxids \
-s 0 \
-t 2 \
-o "$WD"/

exit 0