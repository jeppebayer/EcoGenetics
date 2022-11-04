#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 36:00:00

# Estimation of allele frequency at each site
python people/Jeppe_Bayer/scripts/poolhmm_v1.4.4/pool-hmm.py \
-f people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607 \
-n 100 \
-P 8 \
-S \
-t 0.005

exit 0