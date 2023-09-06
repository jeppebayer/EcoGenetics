#!/bin/bash

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate R
fi

isodir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/02_data_preparation/Isotoma_sp
lepdir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/02_data_preparation/Lepidocyrtus_lignorum

for file in "$isodir"/*/post_filter_stats/*_filtered.coverage; do
    awk -v b="\t" '{print ($1 b $7)}' "$file" > "$(dirname "$file")"/contig.coverage
    Rscript /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/covplot.R "$(dirname "$file")"/contig.coverage
done

for file in "$lepdir"/*/post_filter_stats/*_filtered.coverage; do
    awk -v b="\t" '{print ($1 b $7)}' "$file" > "$(dirname "$file")"/contig.coverage
    Rscript /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/covplot.R "$(dirname "$file")"/contig.coverage
done

exit 0