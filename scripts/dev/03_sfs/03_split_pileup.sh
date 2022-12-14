#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

lines=$1
regionfile=$2

region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$regionfile")

# Function to split pileup file by scaffold
split()
{
    grep "$region" \
    /faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F/Ocin_NYS-F.pileup \
    > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_"$num".mpileup
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

exit 0