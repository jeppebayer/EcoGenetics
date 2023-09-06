#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

export TMPDIR=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/freebayes
mkdir -p "$TMPDIR"

regionfile=$1

region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$regionfile")

# vcf()
# {
#     freebayes \
#     -p 100 \
#     -f /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna \
#     -r "$region" \
#     -v /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".vcf \
#     --pooled-discrete \
#     /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".bam
# }

vcf()
{
    freebayes-parallel \
    <(fasta_generate_regions.py /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna --chunks 100) 5 \
    -p 100 \
    -f /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna \
    -r "$region" \
    -v /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".vcf \
    --pooled-discrete \
    /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_"$num".bam
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