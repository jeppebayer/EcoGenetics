#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
data=$6 # Path to data location in WD
n=$7 # Number of parts
script_path=$8 # Path to script location

# Activates pool-hmm compatible conda environment
source /home/"$USER"/.bashrc
source activate poolhmm

# Funtion for the creation of SFS parts
sfs()
{
    python "$script_path"/modules/poolhmm_v1.4.4/pool-hmm.py \
    -f "$WD"/temp/"$(basename "$sample")"_intergenic_part"$num" \
    -n 100 \
    -P "$cpus" \
    -o \
    -t 0.005 \
    > "$WD"/temp/"$(basename "$sample")"_intergenic_part"$num".spectrum
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
rm -f "$WD"/temp/"$(basename "$sample")"_intergenic_part"$num"

# for part in {0..19..1}; do
#     length=${#part}
#     if [ "$length" -lt 2 ]; then
#         num="0$part"
#         python /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/03_initial_analysis_files/modules/poolhmm_v1.4.4/pool-hmm.py \
#         -f "$WD"/temp/"$(basename "$sample")"_intergenic_part"$num" \
#         -n 100 \
#         -P 10 \
#         -o \
#         -t 0.005 \
#         > "$WD"/temp/"$(basename "$sample")"_intergenic_part"$num".spectrum
#     else
#         num="$part"
#         python /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/03_initial_analysis_files/modules/poolhmm_v1.4.4/pool-hmm.py \
#         -f "$WD"/temp/"$(basename "$sample")"_intergenic_part"$num" \
#         -n 100 \
#         -P 10 \
#         -o \
#         -t 0.005 \
#         > "$WD"/temp/"$(basename "$sample")"_intergenic_part"$num".spectrum
#     fi
# done

exit 0