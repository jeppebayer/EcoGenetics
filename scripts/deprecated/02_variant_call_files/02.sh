#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 36:00:00

# Using pool-hmm script to calculate folded SFS
python people/Jeppe_Bayer/scripts/poolhmm_v1.4.4/pool-hmm.py \
-f people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642 \
-n 100 \
-P 8 \
-o -c 20 -t 0.005

exit 0