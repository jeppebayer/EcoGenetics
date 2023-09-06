#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

target="$1"
WD="$2"
currentuser="$3"

if [ "$USER" == "$currentuser" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly
fi

jellyfish histo \
-t 32 \
-h 100000 \
-o "$WD"/kmer_histogram_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$target")" \
"$WD"/kmer_count_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$target")"

exit 0