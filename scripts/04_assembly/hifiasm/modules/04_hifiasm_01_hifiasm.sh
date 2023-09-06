#!/bin/bash
#SBATCH --account EcoGenetics

targetfile="$1"
WD="$2"
trim="$3"
sim_thres="$4"
speciesabbr="$5"
currentuser="$6"

if [ "$USER" == "$currentuser" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly
fi

if [ -z "$trim" ]; then
    hifiasm \
    -o "$WD"/"$speciesabbr".asm \
    -t 32 \
    -s "$sim_thres" \
    -l 3 \
    "$targetfile"
else
    hifiasm \
    -o "$WD"/"$speciesabbr".trim.asm \
    -t 32 \
    -z "$trim" \
    -s "$sim_thres" \
    -l 3 \
    "$targetfile"
fi

exit 0