#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

target=$1
WD=$2

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

jellyfish histo \
-t 32 \
-o "$WD"/jellyfish_"$(basename "$target")"_histo_"${SLURM_ARRAY_TASK_ID}" \
"$WD"/jellyfish_"$(basename "$target")"_"${SLURM_ARRAY_TASK_ID}"

exit 0