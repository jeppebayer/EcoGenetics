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

split_fa \
"$ref" \
> "$WD"/"$firstletter""$lowerletters".asm.split.fasta

exit 0