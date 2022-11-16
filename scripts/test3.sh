#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

adjustment=0.1

# number=$(awk -v adjustment=$adjustment 'BEGIN { printf "%.0f\n", ( 720 * adjustment) }')
# number=$(awk -v adjustment=$adjustment 'BEGIN { print int( 720 * adjustment + 120) }')
# $((2/filesize))
# echo "$number"

# export DISPLAY=:0

# echo "$PWD"
# echo "$WD"

echo "$(readlink -f "$1")"

exit 0