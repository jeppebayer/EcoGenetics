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

Rscript "$script_path"/modules/genomescope/genomescope.R \
"$WD"/kmer_histogram_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$targetfile")" \
"${SLURM_ARRAY_TASK_ID}" \
"$read_length" \
"$WD"/genomescope_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$targetname")" \
"$kmer_max"

exit 0