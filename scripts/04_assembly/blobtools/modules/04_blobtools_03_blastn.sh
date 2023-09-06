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

# Outfmt: 6 = tabular, additional formatting, std = standard
# Resulting format:
# qseqid (Query Seq-ID), staxids (Unique Subject Taxonomy IDs, separated by a ';'), bitscore (Bit Score), pident (Percentage of identical matches), length (Alignment lenght), mismatch (Number of mismatches), gapopen (Number of gap openings), qstart (Start of alignment in query), qend (End of alignment in query), sstart (Start of alignment in subject), send (End of alignment in subject), evalue (Expect value)
blastn \
-task megablast \
-query "$asmfile" \
-db /faststorage/project/EcoGenetics/BACKUP/database/blastdb/nt \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 10 \
-max_hsps 1 \
-num_threads 16 \
-evalue 1e-25 \
-out "$WD"/"$speciesabbr"_ncbimegablast.out

exit 0