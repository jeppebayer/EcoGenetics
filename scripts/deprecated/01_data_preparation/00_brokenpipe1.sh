#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 48:00:00

# Aligns sample reads to indexed reference genome, pipes into a name sorting the alignment making 
# fixmate possible, then coordinate sorts the alignment and marks duplicates

# Align sample to reference genome and pipe output
bwa mem -t 8 \
BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic \
BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_1.fq.gz BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_2.fq.gz | \

# Sort alignment with regard to QNAME and write output
samtools sort -@ 7 -n -O BAM \
-T people/Jeppe_Bayer/steps/temp2/ \
- | \

# Add fixmate tag to alignment and pipe output. Can only be done on name sorted alignment
samtools fixmate -@ 7 -m -O BAM \
- - | \

# Position sort alignment and pipe output
samtools sort -@ 7 -O BAM \
-T people/Jeppe_Bayer/steps/temp2/ \
- | \

# Mark duplicates and save output as .coordinatesort.bam. Also outputs some realted statistics
samtools markdup -@ 7 -s \
-d 100 \
-f people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.markdupstats \
-T people/Jeppe_Bayer/steps/temp2/ \
- \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/no_intermediate/E100050923_L01_ANIcnqkR228563-642.markdup.bam

exit 0