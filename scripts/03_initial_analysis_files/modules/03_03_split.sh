#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

# Number of files to split to
n=20

# Splits file into n parts without splitting lines. Each part gets samplename + number (00-99) + extension (.part)

split -n l/$n -d --additional-suffix .part input /path/to/temp/"$(basename "$sample")"

exit 0