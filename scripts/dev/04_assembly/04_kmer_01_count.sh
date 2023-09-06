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

if [[ "$target" == *.gz ]]; then
    jellyfish count \
    -t 32 \
    -C \
    -m "${SLURM_ARRAY_TASK_ID}" \
    -s 8G \
    -o "$WD"/jellyfish_"$(basename "$target")"_"${SLURM_ARRAY_TASK_ID}" \
    <(zcat "$target")
else
    jellyfish count \
    -t 32 \
    -C \
    -m "${SLURM_ARRAY_TASK_ID}" \
    -s 8G \
    -o "$WD"/jellyfish_"$(basename "$target")"_"${SLURM_ARRAY_TASK_ID}" \
    "$target"
fi

exit 0