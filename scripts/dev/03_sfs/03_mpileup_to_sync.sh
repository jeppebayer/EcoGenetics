#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

lines=$1
regionfile=$2

region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$regionfile")

# Function to make sync file from mpileup
create_sync()
{
    perl /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/popoolation2_v1.201/mpileup2sync.pl \
    --input /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_filtered_"$num".mpileup \
    --fastq-type sanger \
    --min-qual 0 \
    --output /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_"$num".sync
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

rm -f /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_filtered_"$num".mpileup

exit 0