#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

lines=$1 # Number of unique IDs
qname=$2 # File containing IDs
S=$3 # Path to sample .mpileup file
A=$4 # Number of alleles in sample
namebase=$5 # Current name base for sample
temp=$6 # Directory for temporary files
sampledir=$7 # Sample specific directory within data directory
keep_temp=$8 # Flag for whether or not temporary files should be kept

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary

fi

# Function to make sync file from mpileup
create_sync()
{
    perl /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/popoolation2_v1.201/mpileup2sync.pl \
    --input "$temp"/"$namebase"_filtered_"$num".mpileup \
    --fastq-type sanger \
    --min-qual 1 \
    --output "$temp"/"$namebase"_"$num".sync
}

# Adjusts naming according to file number
id=${SLURM_ARRAY_TASK_ID}
idlength=${#id}
lineslength=${#lines}

if [ $((idlength)) -lt $((lineslength)) ]; then

    dif=$((lineslength - idlength))
    num=""

    for (( i=1;  i<="$dif"; i++ )); do

        num="0$num"

    done

    num="$num$id"
    create_sync

else

    num="$id"
    create_sync
    
fi

if [ "$keep_temp" == "N" ]; then
    rm -f "$temp"/"$namebase"_filtered_"$num".mpileup
fi 


exit 0