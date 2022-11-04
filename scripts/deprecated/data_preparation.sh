#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 48:00:00

bwa mem -t 8 \
BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic \
BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_1.fq.gz BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_2.fq.gz | \
samtools collate -@ 7 \
people/Jeppe_Bayer/steps/temp/ | \
samtools fixmate -@ 7 -m -O BAM | \
samtools sort -@ 7 -O BAM \
-T people/Jeppe_Bayer/steps/temp/ | \
samtools markdup -@ 7 -s \
-f people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_sort/E100050923_L01_ANIcnqkR228563-642.sorted.markdupstats \
-T people/Jeppe_Bayer/steps/temp/ \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam

samtools faidx \
BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna \
> BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna.fai

samtools index -@ 7 -b \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam.bai

samtools flagstat -@ 7 \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_sort/E100050923_L01_ANIcnqkR228563-642.sorted.flagstat

samtools idxstats -@ 7 \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_sort/E100050923_L01_ANIcnqkR228563-642.sorted.idxstats

samtools coverage \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_sort/E100050923_L01_ANIcnqkR228563-642.sorted.coverage

samtools view -b -@ 7 \
-F 0x4 -f 0x2 -q 20 \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam

samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.sorted.bam \
> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/comparison_of_reads.filter && \
samtools view -@ 7 -c \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.filtered.bam \
>> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/comparison_of_reads.filter && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/comparison_of_reads.filter \
>> people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_filter/comparison_of_reads.filter