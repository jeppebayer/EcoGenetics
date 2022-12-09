#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

lines=$1
regionfile=$2

region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$regionfile")

split()
{
    samtools view \
    -@ 9 \
    -b \
    -o /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".bam \
    /faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F/Ocin_NYS-F_filtered.bam \
    "$region"
}
index()
{
    samtools index \
    -@ 9 \
    -b \
    /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".bam \
    > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".bam.bai
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
    index

else

    num="$id"
    split
    index
    
fi

exit 0