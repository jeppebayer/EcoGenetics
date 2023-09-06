#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

targetfile="$1"
WD="$2"
script_path="$3"
kmer_max="$4"
read_length="$5"
currentuser="$6"
targetname="${targetfile%.*}"

if [ "$USER" == "$currentuser" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly
fi

kmerdir="$WD/genomescope2_${SLURM_ARRAY_TASK_ID}_$(basename "$targetname")"
[ -d "$kmerdir " ] || mkdir -m 775 "$kmerdir"

genomescope2 \
-i "$WD"/kmer_histogram_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$targetfile")" \
-p 2 \
-o "$kmerdir" \
-k "${SLURM_ARRAY_TASK_ID}" \
-n "$(basename "$targetname")"_"${SLURM_ARRAY_TASK_ID}" \
-m "$kmer_max"

exit 0