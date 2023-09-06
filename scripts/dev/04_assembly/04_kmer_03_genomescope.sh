#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

target=$1
WD=$2
coverage_max=$3

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

Rscript /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/genomescope/genomescope.R \
"$WD"/jellyfish_"$(basename "$target")"_histo_"${SLURM_ARRAY_TASK_ID}" \
"${SLURM_ARRAY_TASK_ID}" \
100 \
"$WD"/genomescope_"$(basename "$target")"_"${SLURM_ARRAY_TASK_ID}" \
"$coverage_max"

exit 0