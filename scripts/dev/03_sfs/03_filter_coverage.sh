#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

lines=$1
regionfile=$2

region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$regionfile")

# Function to filter mpileup for reads with coverage <= 200
filter()
{
    awk \
    '$4 >= 200' \
    /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_"$num".mpileup \
    > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_filtered_"$num".mpileup
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
    filter

else

    num="$id"
    filter
    
fi

rm -f /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_"$num".mpileup

exit 0