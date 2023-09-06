#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 15G
#SBATCH --cpus-per-task 32
#SBATCH --time 10:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/hifiasm-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

target="$(readlink -f "$1")"
WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly"
name=$(basename "$target")
firstletter=${name:0:1}
nextletters="${name:1:3}"
lowerletters="${nextletters,,}"

# hifiasm \
# -o "$WD"/"$(basename "$(dirname "$target")")"/"$firstletter""$lowerletters".asm \
# -t 32 \
# -s 0.1 \
# -l 3 \
# "$target"

# Trims 20 bp of each end of each contig
hifiasm \
-o "$WD"/"$(basename "$(dirname "$target")")"/"$firstletter""$lowerletters".trim.asm \
-t 32 \
-z 20 \
-s 0.1 \
-l 3 \
"$target"

exit 0