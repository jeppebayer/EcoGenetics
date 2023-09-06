#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    # source activate R
    source activate ecogen_primary
fi

for mapq in /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/02_dev/*/*/*.mapq; do
    if [ -e "$mapq" ]; then
        # Rscript /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/02/02_plotmapq.r "$mapq"
        python /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/02/02_plotmapq.py "$mapq"
    fi
done

exit 0