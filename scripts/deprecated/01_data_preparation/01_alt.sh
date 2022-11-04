#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 36:00:00

AdapterRemoval \
--threads 8 \
--file1 BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_1.fq.gz \
--file2 BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642_2.fq.gz \
--adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
--adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
--minquality 0 \
--minlength 0 \
--basename people/Jeppe_Bayer/data/trimmed_output_paired \
--collapse