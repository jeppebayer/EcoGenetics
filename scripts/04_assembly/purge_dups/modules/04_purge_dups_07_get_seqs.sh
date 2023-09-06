#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

asmfile="$1"
WD="$2"
speciesabbr="$3"
script_path="$4"
currentuser="$5"

if [ "$USER" == "$currentuser" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly
fi

get_seqs \
"$WD"/dups.bed \
"$asmfile"

# Need to make a coverage histogram for the purged file

# Makes coverage histogram
python "$script_path"/modules/hist_plot.py \
-c "$WD"/cutoffs \
"$WD"/PB.stat \
"$WD"/"$speciesabbr"_hist.png

# Creates files with all sequence lengths sorted longest to shortest
bioawk \
-v OFS='\t' \
-c fastx \
'{ print $name, length($seq) }' \
"$WD"/purged.fa \
| sort -k1,1 -k2,2nr \
> "$WD"/sequence_lengths

# Calculates total assembly length
total=$(awk '{Total=Total+$2} END{print Total}' "$WD"/sequence_lengths)
awk '{Total=Total+$2} END{print "Total assembly length is: " Total}' "$WD"/sequence_lengths > "$WD"/asm_lenght_"$total"

exit 0