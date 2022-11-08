#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

# Directory containing scripts (Do NOT end with '/')
scripts="people/Jeppe_Bayer/scripts"

# Species specific reference genome (in FASTA format)
RG="BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna"

# Species specific sample directory (Do NOT end with '/')
SD="BACKUP/population_genetics/collembola/Orchesella_cincta"

# Working directory (Do NOT end with '/')
WD="people/Jeppe_Bayer/steps"

usage()
{
cat << EOF

Usage: 02_00_init_data_prep.sh [-r|--reference] <reference_genome> [-s|--species] <species_directory> [-d|--directory] <working_directory> [-a|--algorithm] <algorithm> [-h|--help]

This script is used for initializing the standardized data preparation procedure for sequence data

PARAMETERS:
    -r | --reference    Species specific reference genome, abosolute path (reference genome in FASTA format)
    -s | --species      Species specific sample directory, abosolute path (Do NOT end with '/')
    -d | --directory    Working directory, abosolute path (Do NOT end with '/')

OPTIONS:
    -a | --algorithm    Choice of algorithm to be used during alignment. mem (>70MB, contemporary samples)[default] or aln (<70MB, historic samples)
    -h | --help         Show this message

EOF
}

RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna"

for file in "$(dirname "$RG")"/*.gff; do
    gff=$file
done

echo "$gff"

exit 1

# path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# path=$(dirname "$path")

# # AdapterRemoval
# jid1=$(sbatch --parsable "$path"/test2.sh)

# Aligning to reference
# sbatch --parsable --dependency=aftany:"$jid1" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/test3.sh

# exit 0