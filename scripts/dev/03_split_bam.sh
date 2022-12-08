#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

regionfile=$1

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
length=${#id}

if [ $(($length)) -lt 2 ]; then
    num="0$id"
    split
    index
else
    num="$id"
    split
    index
fi

exit 0