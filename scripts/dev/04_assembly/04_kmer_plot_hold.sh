#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

target=$1
WD=$2

if [ "$USER" == "jepe" ]; then

    source /home/"$USER"/.bashrc
    source activate R

fi

Rscript \
/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/04_assembly/kmer_plot.r \
"$WD"/jellyfish_"$(basename "$target")"_histo_"${SLURM_ARRAY_TASK_ID}"

exit 0