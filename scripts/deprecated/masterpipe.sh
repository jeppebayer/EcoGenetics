#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH --cpus-per-task 8
#SBATCH --time 48:00:00

bwa mem -t 8 \
reference_genome_prefix \
sample_data1.fq.gz sample_data2.fq.gz | \
samtools view -@ 7 -b | \
samtools sort -@ 7 -O bam -T input.temp \
> output.bam && \
samtools coverage \
output.bam \
> output.coverage