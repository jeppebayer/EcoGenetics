#!/bin/bash

for file in /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/02_dev/*/*/*_markdup.bam; do
    sbatch \
    --time=01:00:00 \
    --mem=64G \
    --cpus-per-task=8 \
    --output="$(dirname "$file")"/stdout_sbatch/"$(basename "$(dirname "$file")")"_unmapped-%j.out \
    /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/02/02_get_unmapped.sh "$file"
done

exit 0