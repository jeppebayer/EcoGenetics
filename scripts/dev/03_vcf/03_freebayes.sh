#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

export TMPDIR=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/freebayes
mkdir -p "$TMPDIR"

lines=$1
regionfile=$2

region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$regionfile")

vcf()
{
    freebayes-parallel \
    <(fasta_generate_regions.py /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna.fai 10000) 10 \
    -p 100 \
    -f /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna \
    -v /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".vcf \
    --pooled-discrete \
    /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".bam
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
    vcf

else
    
    num="$id"
    vcf

fi


exit 0