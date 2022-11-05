#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

# ------------------------------------------------------------------------
# Created by Jeppe Bayer
# 
# Master's student at Aarhus University
# Department of Genetics, Ecology and Evolution
# Centre for EcoGenetics
# 
# For troubleshooting or questions:
# Email: jeppe.bayer@bio.au.dk
# ------------------------------------------------------------------------

# ----------------- Description ------------------------------------------

# Script for initializing data preparation procedure for sequence data >70MB

# ----------------- Configuration ----------------------------------------

# Directory containing scripts, abosolute path (Do NOT end with '/')
scripts="/home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/02_data_preparation"

# Species specific reference genome, abosolute path (in FASTA format)
RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna"

# Species specific sample directory, abosolute path (Do NOT end with '/')
SD="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta"

# Working directory, abosolute path (Do NOT end with '/')
WD="/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps"

# ----------------- Script Queue -----------------------------------------

# Creates temp directory in working directory if none exist
[[ -d $WD/temp ]] || mkdir -m 764 $WD/temp

# Creates 01_data_preparation directory in working directory if none exist
[[ -d $WD/01_data_preparation ]] || mkdir -m 764 $WD/01_data_preparation

# Creates species directory in 01_data_preparation directory if none exist
[[ -d $WD/01_data_preparation/$(basename $SD) ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename $SD)"

# Loops through all sample folders within species specific sample directory
for sample in "$SD"/*; do

    # Checks if sample folder is empty
    if [ "$(ls -A "$sample")" ]; then

        for file in "$sample"/*.bam; do

            # Checks whether a .bam file already exists within sample folder, indicating samples have already been processed
            if [ ! -e "$file" ]; then
                
                # Creates sample directory in species directory if none exist
                [[ -d $WD/01_data_preparation/$(basename $SD)/"$(basename "$sample")" ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename $SD)"/"$(basename "$sample")"
                
                # Creates pre- and post-filtering directory in sample directory if none exist
                [[ -d $WD/01_data_preparation/$(basename $SD)/"$(basename "$sample")"/pre_filter_stats ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename $SD)"/"$(basename "$sample")"/pre_filter_stats
                [[ -d $WD/01_data_preparation/$(basename $SD)/"$(basename "$sample")"/post_filter_stats ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename $SD)"/"$(basename "$sample")"/post_filter_stats
                
                # AdapterRemoval
                jid1=$(sbatch --parsable "$scripts"/02_01_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                # Aligning to reference
                jid2_1=$(sbatch --parsable --dependency=afterany:"$jid1" "$scripts"/02_02_01_data_prep.sh "$RG" "$SD" "$WD" "$sample")
                jid2_2=$(sbatch --parsable --dependency=afterany:"$jid1" "$scripts"/02_02_02_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                # Merging of alignment files
                jid3=$(sbatch --parsable --dependency=afterany:"$jid2_1":"$jid2_2" "$scripts"/02_03_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                # Marking duplicates
                jid4=$(sbatch --parsable --dependency=afterany:"$jid3" "$scripts"/02_04_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                # Statistics pre-filtering
                jid5_1=$(sbatch --parsable --dependency=afterany:"$jid4" "$scripts"/02_05_01_data_prep.sh "$RG" "$SD" "$WD" "$sample")
                jid5_2=$(sbatch --parsable --dependency=afterany:"$jid4" "$scripts"/02_05_02_data_prep.sh "$RG" "$SD" "$WD" "$sample")
                jid5_3=$(sbatch --parsable --dependency=afterany:"$jid4" "$scripts"/02_05_03_data_prep.sh "$RG" "$SD" "$WD" "$sample")
                jid5_4=$(sbatch --parsable --dependency=afterany:"$jid4" "$scripts"/02_05_04_data_prep.sh "$RG" "$SD" "$WD" "$sample")


                # Removal of duplicates, unmapped reads and low quality mappings
                jid6=$(sbatch --parsable --dependency=afterany:"$jid5_1":"$jid5_2":"$jid5_3":"$jid5_4" "$scripts"/02_06_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                # Statistics post-filtering
                sbatch --dependency=afterany:"$jid6" "$scripts"/02_07_01_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                sbatch --dependency=afterany:"$jid6" "$scripts"/02_07_02_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                sbatch --dependency=afterany:"$jid6" "$scripts"/02_07_03_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                sbatch --dependency=afterany:"$jid6" "$scripts"/02_07_04_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                sbatch --dependency=afterany:"$jid6" "$scripts"/02_07_05_data_prep.sh "$RG" "$SD" "$WD" "$sample"

            else

                srun echo "$sample already contains a .bam file, $file, and is skipped"
            fi

        done

    else

        srun echo "$sample is an empty directory and is skipped"
    fi

done

exit 0