#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --chdir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/blastdb
#SBATCH --mem-per-cpu 12G
#SBATCH --cpus-per-task 5
#SBATCH --time 08:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/blasdb_download-%j.out

wget -q -i URLs.txt

exit 0