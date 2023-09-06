#!/bin/bash

datapath=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/X_evolution/data/individual_data

for original_file in "$datapath"/*/*/*filtered.bam*; do
    rm "$original_file"
done
for original_file in "$datapath"/*/*/*unmapped.bam*; do
    rm "$original_file"
done

for retag_file in "$datapath"/*/*/*retag*; do
    extension="${retag_file#*.}"
    name="${retag_file%*_retag*}"
    mv "$retag_file" "$name"."$extension"
done

exit 0