#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

for bam in /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/02_dev/*/*/*_markdup.bam; do
    for mapq in "$(dirname "$bam")"/*.mapq; do
        if [ ! -e "$mapq" ]; then
            samtools view "$bam" | awk -F '\t' '{print $5}' - > "$(dirname "$bam")"/"$(basename "$(dirname "$bam")")".mapq
        fi
    done
done

exit 0