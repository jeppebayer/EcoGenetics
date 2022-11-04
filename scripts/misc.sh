#!/bin/bash

# Path to reference
reference=/home/jepe/EcoGenetics/BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1/
# Reference file suffix
rsuffix=.fna

# Path to folder containing sample
sample=/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F/
# Sample file suffix
ssuffix=.fq.gz

# Work folder
wfolder=/home/jepe/EcoGenetics/people/Jeppe_Bayer/

test -d people/Jeppe_Bayer/data || mkdir -m 0764 people/Jeppe_Bayer/data

exit 0