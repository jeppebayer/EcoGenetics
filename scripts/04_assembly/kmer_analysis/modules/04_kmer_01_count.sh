#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

target="$1"
WD="$2"
currentuser="$3"
canonical="$4"

if [ "$USER" == "$currentuser" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly
fi

if [ -z "$canonical" ]; then
    if [[ "$target" == *.gz ]]; then
        jellyfish count \
        -t 32 \
        -m "${SLURM_ARRAY_TASK_ID}" \
        -s 8G \
        -o "$WD"/kmer_count_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$target")" \
        <(zcat "$target")
    else
        jellyfish count \
        -t 32 \
        -m "${SLURM_ARRAY_TASK_ID}" \
        -s 8G \
        -o "$WD"/kmer_count_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$target")" \
        "$target"
    fi
else
    if [[ "$target" == *.gz ]]; then
        jellyfish count \
        -t 32 \
        -C \
        -m "${SLURM_ARRAY_TASK_ID}" \
        -s 8G \
        -o "$WD"/kmer_count_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$target")" \
        <(zcat "$target")
    else
        jellyfish count \
        -t 32 \
        -C \
        -m "${SLURM_ARRAY_TASK_ID}" \
        -s 8G \
        -o "$WD"/kmer_count_"${SLURM_ARRAY_TASK_ID}"_"$(basename "$target")" \
        "$target"
    fi
fi

exit 0