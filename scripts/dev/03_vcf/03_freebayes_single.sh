#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --time=240
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=10
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/out/freebayes-test13-%j.out

export TMPDIR=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/freebayes
mkdir -p "$TMPDIR"

regionfile=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/qname.txt

region=$(sed -n "1p" "$regionfile")

freebayes-parallel \
<(fasta_generate_regions.py /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna.fai 10000) 10 \
-p 100 \
-f /faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna \
-v /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/a_test13.vcf \
-n 4 \
--pooled-discrete \
/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_01.bam