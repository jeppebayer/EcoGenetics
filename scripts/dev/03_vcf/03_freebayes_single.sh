#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --chdir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp
#SBATCH --time=1200
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/a_test19-%j.vcf

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

export TMPDIR=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/freebayes
[ -d "$TMPDIR" ] || mkdir -m 775 "$TMPDIR"

regionfile=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/qname.txt

region=$(sed -n "1p" "$regionfile")

# freebayes-parallel \
# <(fasta_generate_regions.py /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna.fai 10000) 2 \
# -p 100 \
# -f /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna \
# -v /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/a_test15.vcf \
# -n 3 \
# --pooled-discrete \
# /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_01.bam

# freebayes \
# -p 100 \
# -f /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna \
# -v /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/a_test17.vcf \
# -n 5 \
# -r "$region" \
# --pooled-discrete \
# /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_01.bam

freebayes-parallel \
<(fasta_generate_regions.py /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna.fai --chunks 100) 8 \
-p 100 \
-f /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna \
-r "$region" \
--pooled-discrete \
/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_01.bam

exit 0
