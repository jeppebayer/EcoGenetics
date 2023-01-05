#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

target="$1"
ref="$2"
WD="$3"
name=$(basename "$target")
firstletter=${name:0:1}
nextletters="${name:1:3}"
lowerletters="${nextletters,,}"

get_seqs \
"$WD"/dups.bed \
"$ref"

# Makes coverage histogram
python /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/04_assembly/hist_plot.py \
-c "$WD"/cutoffs \
"$WD"/PB.stat \
"$WD"/"$firstletter""$lowerletters"_hist.png

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