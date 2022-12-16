#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 00:40:00
#SBATCH --output=assembly-%j.out

temp="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp"
scripts="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/04_assembly"
WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly"
out="$WD"/out

[ -d "$temp" ] || mkdir -m 775 "$temp"
[ -d "$WD" ] || mkdir -m 775 "$WD"
[ -d "$out" ] || mkdir -m 775 "$out"

jid1=$(sbatch \
    --parsable \
    --time=120 \
    --mem-per-cpu=16G \
    --cpus-per-task=2 \
    --output="$out"/index-%j.out \
    "$scripts"/04_index.sh)