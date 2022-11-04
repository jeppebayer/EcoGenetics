#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 36:00:00

# Align sample to reference genome
bwa mem -t 8 \
BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic \
BACKUP/population_genetics/collembola/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607_1.fq.gz BACKUP/population_genetics/collembola/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607_2.fq.gz \
| \

# Sort alignment with regard to QNAME
samtools sort -@ 7 -n -O BAM \
-T people/Jeppe_Bayer/steps/temp/ \
- | \

# Add fixmate tag to alignment and pipe output. Can only be done on name sorted alignment
samtools fixmate -@ 7 -m -O BAM \
- - | \

# Position sort alignment and pipe output
samtools sort -@ 7 -O BAM \
-T people/Jeppe_Bayer/steps/temp/ \
- | \

# Creates bam file of name sorted reads
tee \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.positionsort.bam \
| \

# Mark duplicates and save output as .coordinatesort.bam. Also outputs some realted statistics
samtools markdup -@ 7 -s \
-f people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/QC_sort/E100049659_L01_ANIcnqkR228559-607.markdup.markdupstats \
-T people/Jeppe_Bayer/steps/temp/ \
- \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.markdup.bam

exit 0