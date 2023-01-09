#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --chdir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/blastdb
#SBATCH --mem-per-cpu 12G
#SBATCH --cpus-per-task 2
#SBATCH --time 01:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/extract-%j.out

# for file in ./*.tar.gz; do
    
#     tar -zxf "$file" -C .

# done

for file in ./*.tar.gz; do
    
    rm -f "$file"

done


exit 0