#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 6
#SBATCH --time 20:00:00

samtools faidx \
BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna \
> BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna.fai && \
samtools index -@ 7 -b \
people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/E100050923_L01_ANIcnqkR228563-642.bam \
> people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/E100050923_L01_ANIcnqkR228563-642.bam.bai && \
platypus callVariants \
--bamFile=people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/E100050923_L01_ANIcnqkR228563-642.bam \
--refFile=BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna \
--output=people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/filtered/E100050923_L01_ANIcnqkR228563-642.vcf