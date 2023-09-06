#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

bam="$1"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ 7 \
-f 4 \
-o "$(dirname "$bam")"/"$(basename "$(dirname "$bam")")"_unmapped.bam \
"$bam" && \

# Creates bai index for filtered alignment
samtools index -@ 7 -b \
"$(dirname "$bam")"/"$(basename "$(dirname "$bam")")"_unmapped.bam \
> "$(dirname "$bam")"/"$(basename "$(dirname "$bam")")"_unmapped.bam.bai

exit 0