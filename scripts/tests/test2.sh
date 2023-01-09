#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --chdir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/uniprot
#SBATCH --mem-per-cpu 12G
#SBATCH --cpus-per-task 2
#SBATCH --time 10:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/unzip_and_database-%j.out

# wget -q https://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz

gzip -d uniprot_trembl.fasta.gz

diamond makedb --in uniprot_trembl.fasta -d uniprot_trembl

exit 0