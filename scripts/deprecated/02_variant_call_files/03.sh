#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH --cpus-per-task 8
#SBATCH --time 36:00:00

# Estimation of allele frequency at each site
python people/Jeppe_Bayer/scripts/poolhmm_v1.4.4/pool-hmm.py \
-f people/Jeppe_Bayer/data/Ocin_NYS-F_6 \
-n 100 \
-P 8 \
-o \
-t 0.005

exit 0