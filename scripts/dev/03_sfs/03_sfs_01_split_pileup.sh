#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

export LC_ALL=C

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

region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$qname")

# Function to split pileup file by scaffold
split()
{
    rg \
    -j 6 \
    -F \
    "$region" \
    "$S" \
    > "$temp"/"$namebase"_"$num".mpileup
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
    split

else

    num="$id"
    split
    
fi

if [ "$keep_temp" == "N" ]; then
    rm -f "$qname"
fi 

exit 0