#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

targetfile="$1"
currentuser="$2"
fasta_name="${targetfile%.*}"

if [ "$USER" == "$currentuser" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

awk '/^S/{print ">"$2"\n"$3}' "$targetfile" | fold > "$fasta_name".fasta

# Creates files with all sequence lengths sorted longest to shortest
bioawk \
-v OFS='\t' \
-c fastx \
'{ print $name, length($seq) }' \
"$fasta_name".fasta \
| sort -k1,1 -k2,2nr \
> "$fasta_name"_sequence_lengths

# Calculates total assembly length
total=$(awk '{Total=Total+$2} END{print Total}' "$fasta_name"_sequence_lengths)
# awk '{Total=Total+$2} END{print "Total assembly length is: " Total}' "$fasta_name"_sequence_lengths > "$fasta_name"_asm_lenght_"$total"
echo "Total assembly length is: $total" > "$fasta_name"_asm_lenght_"$total"

exit 0