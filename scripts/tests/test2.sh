#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 1
#SBATCH --time 01:00:00

# awk '$4 >= 200' /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test2.txt

# Removes empty .out files
for file in /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs/out/*; do
    if [ ! -s "$file" ]; then
        rm -f "$file"
    fi
done

exit 0