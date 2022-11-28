#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 30:00:00
#SBATCH --output=/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps/03_initial_analysis_files/03_sfs_full-%j.out

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
data=$6 # Path to data location in WD

# Activates pool-hmm compatible conda environment
source /home/"$USER"/.bashrc
source activate poolhmm

RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna"
sample="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F"
WD="/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps"

# Funtion for the creation of SFS parts
sfs()
{
    python /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/03_initial_analysis_files/modules/poolhmm_v1.4.4/pool-hmm.py \
    -f "$WD"/temp/"$(basename "$sample")"_nonsyn_part"$num" \
    -n 100 \
    -P "$cpus" \
    -o \
    -t 0.005 \
    > "$WD"/temp/"$(basename "$sample")"_nonsyn_part"$num".spectrum
}

# Adjusts naming according to file number
length=${#$SLURM_ARRAY_TASK_ID}

if [ "$length" -lt 2 ]; then
    num="0$SLURM_ARRAY_TASK_ID"
    sfs
else
    num="$SLURM_ARRAY_TASK_ID"
    sfs
fi

# File removal
rm -f "$WD"/temp/"$(basename "$sample")"_nonsyn_part"$num"

exit 0