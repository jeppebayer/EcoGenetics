#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:40:00
#SBATCH --output=Init_Analysis_Files_Initialize-%j.out

# ----------------- Description ------------------------------------------

# Script for initializing the creation of some standard files for further analyses
# Tested and working using EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 03_initialize_init_analysis_files.sh [PARAMETERS] [OPTIONS]

This script is used for initializing the creation of some standard files for further analyses.
Intended to be used in conjunction with 'sbatch', however it also be used in conjunction with 'srun' 
or simply on the frontend as its resource-demands are fairly low.

PARAMETERS (must be assigned):
    -r  FILE            Species specific reference genome, abosolute path (reference genome in FASTA format)
    -s  DIRECTORY       Species specific sample directory, abosolute path (Do NOT end with '/')
    -d  DIRECTORY       Working directory, abosolute path (Do NOT end with '/')

OPTIONS:
    -m  INTEGER         Amount of memory to be used by each CPU. 8 [default]
    -c  INTEGER         Number of CPUs to be used. 8 [default]
    -h                  Show this message

Parameters can also be set directly in the script. 
If you decide to do so, change the variables under CONFIGURATION.

Tested and working using 'EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml'

EOF
}

# ----------------- Configuration ----------------------------------------

# Species specific reference genome, abosolute path (reference genome in FASTA format)
RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna"

# Species specific sample directory, abosolute path (Do NOT end with '/')
SD="/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta"

# Working directory, abosolute path (Do NOT end with '/')
WD="/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps"

# Define memory per cpu in G (must be integer)
memory="8"

# Define number of cpus to be used (must be integer)
cpus="8"

# ----------------- Error messages ---------------------------------------

# Gets path to script location
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

error_r()
{
cat << EOF

ERROR: $OPTARG does not seem to be a .fna file

If unsure of how to proceed run: $(basename "$script_path") -h

EOF
}

error_s()
{
cat << EOF

ERROR: $OPTARG does not seem to be a directory

If unsure of how to proceed run: $(basename "$script_path") -h

EOF
}

error_d()
{
cat << EOF

ERROR: $OPTARG does not seem to be a directory

If unsure of how to proceed run: $(basename "$script_path") -h

EOF
}

error_m()
{
cat << EOF

ERROR: -m must be a whole integer, 8 [default]'

If unsure of how to proceed run: $(basename "$script_path") -h

EOF
}

error_c()
{
cat << EOF

ERROR: -c must be a whole integer, 8 [default]'

If unsure of how to proceed run: $(basename "$script_path") -h

EOF
}

# ----------------- Script Flag Processing -------------------------------

while getopts 'r:s:d:m:c:h' OPTION; do
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
        m)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                memory="$OPTARG"
            else
                error_m
                exit 1
            fi
            ;;
        c)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                cpus="$OPTARG"
            else
                error_c
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
[[ -d "$WD" ]] || mkdir -m 775 "$WD"

# Creates temp directory in working directory if none exist
[[ -d "$WD"/temp ]] || mkdir -m 775 "$WD"/temp

# Creates 01_data_preparation directory in working directory if none exist
[[ -d "$WD"/"$(basename "$script_path")" ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"

# Creates species directory in 01_data_preparation directory if none exist
[[ -d "$WD"/"$(basename "$script_path")"/"$(basename "$SD")" ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"

# Creates log file
touch "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

for sample in "$SD"/*; do

    # Checks if sample folder is empty
    if [ "$(ls -A "$sample")" ]; then

        # Checks whether a .bam file exists within sample folder
        for file in "$sample"/*.bam; do
            
            if [ -e "$file" ]; then

                # Creates sample directory in species directory if none exist
                [[ -d "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")" ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"

                # Create folder for STDOUT files generated sbatch
                [[ -d "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch
                
                echo "$(basename "$sample") is sent to the queue" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                # Make mpileup file
                jid1=$(sbatch \
                            --parsable \
                            --time=1440 \
                            --mem-per-cpu="$memory"G \
                            --cpus-per-task="$cpus" \
                            --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_01-%j.out \
                            "$script_path"/modules/03_01_mpileup.sh "$RG" "$SD" "$WD" "$sample")

                echo "'mpileup' job has been submitted for $(basename "$sample") -- Job ID: $jid1" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                # Make VCF file
                jid2=$(sbatch \
                            --parsable \
                            --time=1440 \
                            --mem-per-cpu="$memory"G \
                            --cpus-per-task="$cpus" \
                            --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_01-%j.out \
                            --dependency=afterany:"$jid1" \
                            "$script_path"/modules/03_02_vcf.sh "$RG" "$SD" "$WD" "$sample")                

                echo "'VCF' job has been submitted for $(basename "$sample") -- Job ID: $jid2" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

            else

                echo "$(basename "$sample") doesn't contain a .bam file and is skipped" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

            fi

        done

    else

        echo "$(basename "$sample") is an empty directory and is skipped" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

    fi

done

echo "Log file can be found here: $WD/$(basename "$script_path")/$(basename "$SD")/log_$(basename "$SD")_$(date +"%Y-%m-%d").txt"

exit 0
