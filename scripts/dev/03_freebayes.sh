#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

regionfile=$1

region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$regionfile")

vcf()
{
    freebayes \
    -p 100 \
    -f /home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna \
    -r "$region" \
    -v /home/jepe/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".vcf \
    --pooled-discrete \
    /home/jepe/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".bam
}

# Adjusts naming according to file number
id=${SLURM_ARRAY_TASK_ID}
length=${#id}

if [ $(($length)) -lt 2 ]; then
    num="0$id"
    vcf
else
    num="$id"
    vcf
fi

exit 0