#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --time=00:30:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=8
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/dev_depth-%j.out

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

bamfile=$(readlink -f "$1")
outpath=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/02_data_preparation/"$(basename "$(dirname "$(dirname "$bamfile")")")"/"$(basename "$(dirname "$bamfile")")"/post_filter_stats

samtools depth \
-@ 7 \
-a \
-o "$outpath/$(basename "$(dirname "$bamfile")")_filtered.depth" \
"$bamfile"

exit 0