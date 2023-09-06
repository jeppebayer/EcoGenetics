#!/bin/bash

gtf_file=/faststorage/project/EcoGenetics/people/EliseL/Functional_neutral_GD/Populations/Outputs/Orchesella_cincta/Orc-cinc_functionalGD_particulier.gtf
dist_file=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/functional_dist.txt

gene_id=($(awk \
    -F "\t" \
    '{print $9}' \
    "$gtf_file" \
| awk \
    -F ";" \
    '{print $1}' \
| awk \
    '{print $2}' \
| uniq))

echo -n "" > "$dist_file"

for gene in "${gene_id[@]}"; do
    length=$(grep "$gene" "$gtf_file" \
    | awk \
        -F "\t" \
        'BEGIN {OFS="\t"}
        {Total += ($5 - $4)}
        END {print Total}')
    echo -e "$gene\t$length" >> "$dist_file"
done

exit 0