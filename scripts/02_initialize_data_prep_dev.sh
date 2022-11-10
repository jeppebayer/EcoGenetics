#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:40:00

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

# Script for initializing data preparation procedure for sequence data
# Tested and working using EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 02_initialize_data_prep.sh [-r FILE] [-s DIRECTORY] [-d DIRECTORY] [-a ALGORITHM] [-h]

This script is used for initializing the standardized data preparation procedure for sequence data.
Intended to be used in conjunction with 'sbatch', however it also be used in conjunction with 'srun' 
or simply on the frontend as its resource-demands are fairly low.

The 'Mads loop condition': If the specified species sample directory is within the 'museomics' directory, 
the script will automatically try to detect whether the sample directory is contemporary or historical 
and automatically applies the corresponding algorithm. This means that if your samples are in the 
museomics directory you do not need to specify an algorithm.

PARAMETERS:
    -r  FILE            Species specific reference genome, abosolute path (reference genome in FASTA format)
    -s  DIRECTORY       Species specific sample directory, abosolute path (Do NOT end with '/')
    -d  DIRECTORY       Working directory, abosolute path (Do NOT end with '/')

OPTIONS:
    -a  ALGORITHM       Choice of algorithm to be used during alignment. 'mem' (>70MB, contemporary samples)[default] or 'aln' (<70MB, historic samples)
    -h                  Show this message

Parameters can also be set directly in the script. 
If you decide to do so, change the variables under CONFIGURATION.

Tested and working using 'EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml'

EOF
}

# ----------------- Configuration ----------------------------------------

# Species specific reference genome, abosolute path (reference genome in FASTA format)
RG=

# Species specific sample directory, abosolute path (Do NOT end with '/')
SD=

# Working directory, abosolute path (Do NOT end with '/')
WD=

# Algorithm to use during alignment: mem (>70MB, contemporary samples) [default] or aln (<70MB, historic samples)
algo="mem"

# ----------------- Error messages ---------------------------------------

error_r()
{
cat << EOF

ERROR: $OPTARG does not seem to be a .fna file

If unsure of how to proceed run: 02_initialize_data_prep.sh -h

EOF
}

error_s()
{
cat << EOF

ERROR: $OPTARG does not seem to be a directory

If unsure of how to proceed run: 02_initialize_data_prep.sh -h

EOF
}

error_d()
{
cat << EOF

ERROR: $OPTARG does not seem to be a directory

If unsure of how to proceed run: 02_initialize_data_prep.sh -h

EOF
}

error_a()
{
cat << EOF

ERROR: -a must be either 'mem' [default] or 'aln'

If unsure of how to proceed run: 02_initialize_data_prep.sh -h

EOF
}

# ----------------- Script Flag Processing -------------------------------

while getopts 'r:s:d:a:h' OPTION; do
    case "$OPTION" in
        r)
            if [ -f "$OPTARG" ] && [[ $OPTARG == *.fna ]]; then
                RG="$OPTARG"
            else
                error_r
                exit 1
            fi
            ;;
        s)
            if [ -d "$OPTARG" ]; then
                SD="$OPTARG"
            else
                error_s
                exit 1
            fi
            ;;
        d)
            if [ -d "$OPTARG" ]; then
                WD="$OPTARG"
            else
                error_d
                exit 1
            fi
            ;;
        a)
            if [[ $OPTARG == "mem" ]] || [[ $OPTARG == "aln" ]]; then
                algo="$OPTARG"
            else
                error_a
                exit 1
            fi
            ;;
        h)
            usage
            exit 1
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done

if [[ -z $RG ]] || [[ -z $SD ]] || [[ -z $WD ]]; then
    usage
    exit 1
fi

# ----------------- Script Queue -----------------------------------------

# Gets path to script location
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
    script_path=$(dirname "$script_path")
# If run on the frontend:
else
    script_path=$(realpath "$0")
    script_path=$(dirname "$script_path")
fi

# Creates working directory if it doesn't exist
[[ -d "$WD" ]] || mkdir -m 764 "$WD" # || mkdir -m 764 "$PWD"/"$WD"

# Creates temp directory in working directory if none exist
[[ -d "$WD"/temp ]] || mkdir -m 764 "$WD"/temp

# Creates 01_data_preparation directory in working directory if none exist
[[ -d "$WD"/01_data_preparation ]] || mkdir -m 764 "$WD"/01_data_preparation

# Creates species directory in 01_data_preparation directory if none exist
[[ -d "$WD"/01_data_preparation/"$(basename "$SD")" ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename "$SD")"

# Creates log file
touch "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

# Loops through all sample folders within species specific sample directory
for sample in "$SD"/*; do

    # Checks if sample folder is empty
    if [ "$(ls -A "$sample")" ]; then
        
        # Checks whether a .bam file already exists within sample folder, indicating samples have already been processed        
        for file in "$sample"/*.bam; do
            if [ ! -e "$file" ]; then
                
                # Checks sample file size
                filesize=
                for file in "$sample"/*.fna; do
                    if [ ! "$filesize" ]; then
                        filesize=$(wc -c < "$file")
                    fi
                done

                # Change requested time on nodes depending on filesize comparative to R1 Ocin_NYS-F
                adjustment=$(awk -v filesize=$filesize 'BEGIN { print ( filesize / 93635424798 ) }')

                # Creates sample directory in species directory if none exist
                [[ -d "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")" ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"

                # Creates pre- and post-filtering directory in sample directory if none exist
                [[ -d "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats
                [[ -d "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats

                # Checks if currently working with museomics samples
                if [[ "$SD" == *"museomics"* ]]; then

                    # Checks if sample directory is pre-2000 (historic) or post-2000 (modern) and sets the algrorithm parameter correspondingly
                    if [ "$((${sample: -4}))" -gt 2000 ]; then

                        # Sample is modern
                        algo="mem"

                        # Check the number of .fq.gz files in sample directory (assumed to be indicative of whether sample is single- or paired-end)
                        count=$(find "$sample"/ -maxdepth 1 -type f -name '*.fq.gz' | wc -l)
                        if [ "$count" == 1 ]; then
                            
                            # Single-end
                            
                            # AdapterRemoval
                            jid1=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 720 * adjustment + 120) }') "$script_path"/modules/02_01_single_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                            # Aligning to reference
                            jid2=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 1800 * adjustment + 120) }') --dependency=afterany:"$jid1" "$script_path"/modules/02_02_single_data_prep.sh "$RG" "$SD" "$WD" "$sample" "$algo")

                            # Marking duplicates
                            jid4=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 1800 * adjustment + 120) }')--dependency=afterany:"$jid2" "$script_path"/modules/02_04_data_prep.sh "$RG" "$SD" "$WD" "$sample")
                        
                        else

                            # Paired-end
                            
                            # AdapterRemoval
                            jid1=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 720 * adjustment +120 ) }') "$script_path"/modules/02_01_paired_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                            # Aligning to reference
                            jid2_1=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 1800 * adjustment +120 ) }') --dependency=afterany:"$jid1" "$script_path"/modules/02_02_paired_01_data_prep.sh "$RG" "$SD" "$WD" "$sample" "$algo")
                            jid2_2=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 1800 * adjustment +120 ) }') --dependency=afterany:"$jid1" "$script_path"/modules/02_02_paired_02_data_prep.sh "$RG" "$SD" "$WD" "$sample" "$algo")

                            # Merging of alignment files
                            jid3=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 1200 * adjustment +120 ) }') --dependency=afterany:"$jid2_1":"$jid2_2" "$script_path"/modules/02_03_paired_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                            # Marking duplicates
                            jid4=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 1800 * adjustment +120 ) }') --dependency=afterany:"$jid3" "$script_path"/modules/02_04_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                        fi
                        
                        # Statistics pre-filtering
                        jid5_1=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 60 * adjustment + 60) }') --dependency=afterany:"$jid4" "$script_path"/modules/02_05_01_data_prep.sh "$RG" "$SD" "$WD" "$sample")
                        jid5_2=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 60 * adjustment + 60) }') --dependency=afterany:"$jid4" "$script_path"/modules/02_05_02_data_prep.sh "$RG" "$SD" "$WD" "$sample")
                        jid5_3=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 360 * adjustment + 120 ) }') --dependency=afterany:"$jid4" "$script_path"/modules/02_05_03_data_prep.sh "$RG" "$SD" "$WD" "$sample")
                        jid5_4=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 60 * adjustment + 60) }') --dependency=afterany:"$jid4" "$script_path"/modules/02_05_04_data_prep.sh "$RG" "$SD" "$WD" "$sample")


                        # Removal of duplicates, unmapped reads and low quality mappings
                        jid6=$(sbatch --parsable --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 1440 * adjustment + 120 ) }') --dependency=afterany:"$jid5_1":"$jid5_2":"$jid5_3":"$jid5_4" "$script_path"/modules/02_06_data_prep.sh "$RG" "$SD" "$WD" "$sample")

                        # Statistics post-filtering
                        sbatch --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 60 * adjustment + 60 ) }') --dependency=afterany:"$jid6" "$script_path"/modules/02_07_01_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                        sbatch --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 60 * adjustment + 60 ) }') --dependency=afterany:"$jid6" "$script_path"/modules/02_07_02_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                        sbatch --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 360 * adjustment + 120 ) }') --dependency=afterany:"$jid6" "$script_path"/modules/02_07_03_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                        sbatch --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 360 * adjustment + 120 ) }') --dependency=afterany:"$jid6" "$script_path"/modules/02_07_04_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                        sbatch --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 60 * adjustment + 60 ) }') --dependency=afterany:"$jid6" "$script_path"/modules/02_07_05_data_prep.sh "$RG" "$SD" "$WD" "$sample"
                        sbatch --time=$(awk -v adjustment=$adjustment 'BEGIN { print int( 1200 * adjustment + 120 ) }') --dependency=afterany:"$jid6" "$script_path"/modules/02_07_06_data_prep.sh "$RG" "$SD" "$WD" "$sample"

                        # Log messages
                        if [ "$count" == 1 ]; then
                            echo "$sample has been sent to queue as a single-end modern sample" >> "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        else
                            echo "$sample has been sent to queue as a pair-end modern sample" >> "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        fi

                    else

                        # Sample is historic
                        algo="aln"

                        # Check the number of .fq.gz files in sample directory (assumed to be indicative of whether sample is single- or paired-end)
                        count=$(find "$sample"/ -maxdepth 1 -type f -name '*.fq.gz' | wc -l)
                        if [ "$count" == 1 ]; then
                            
                            # Single-end
                        
                        else

                            # Paired-end

                        fi

                        # Log messages
                        if [ "$count" == 1 ]; then
                        
                            echo "$sample has been sent to queue as a single-end historic sample" >> "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        else
                            echo "$sample has been sent to queue as a pair-end historic sample" >> "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        fi

                    fi

                else

                    # Check the number of .fq.gz files in sample directory (assumed to be indicative of whether sample is single- or paired-end)
                    count=$(find "$sample"/ -maxdepth 1 -type f -name '*.fq.gz' | wc -l)
                    if [ "$count" == 1 ]; then

                        # Single-end
                        
                    else

                        # Paired-end

                    fi

                    # Log messages
                    if [ "$count" == 1 ]; then
                        echo "$sample has been sent to queue as a single-end sample" >> "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                    else
                        echo "$sample has been sent to queue as a pair-end sample" >> "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                    fi

                fi

            else

                echo "$sample already contains a .bam file, $file, and is skipped" >> "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

            fi

        done

    else

        echo "$sample is an empty directory and is skipped" >> "$WD"/01_data_preparation/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

    fi

done

echo "Log file can be found here: $WD/01_data_preparation/$(basename "$SD")/log_$(basename "$SD")_$(date +"%Y-%m-%d").txt"

exit 0

93635424798
5438628585

