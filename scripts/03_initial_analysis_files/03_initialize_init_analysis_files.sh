#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:40:00
#SBATCH --output=Init_Analysis_Files_Initialize-%j.out

# ------------------------------------------------------------------------
# Created by Jeppe Bayer
# 
# Master's student at Aarhus University
# Department of Genetics, Ecology and Evolution
# Centre for EcoGenetics
# 
# For troubleshooting or questions:
# Email: jeppe.bayer@bio.au.dk

# ----------------- Description ------------------------------------------

# Script for creating common initial data files to be used for further analyses
# Tested and working using:
# EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 03_initialize_init_analysis_files.sh [PARAMETERS] [OPTIONS]

This script is used for the creation of a series of common initial data files
to be used for further analyses of sequence data at the Center for EcoGenetics.
Can be run on the frontend, as its resource-demands are fairly, or run in 
conjunction with 'sbatch'or 'srun'.

PARAMETERS (must be assigned):
    -s  DIRECTORY       Species specific sample directory

OPTIONS:
    -d  DIRECTORY       Working directory. If not assigned, will use current 
                        working directory [default]
    -m  INTEGER         Amount of memory to be used by each CPU. 8 [default]
    -c  INTEGER         Number of CPUs to be used. 8 [default]
    -f                  Force run even if target sample directory contains relevant files
    -h                  Show this message

Parameters can also be set directly in the script. 
If you decide to do so, change the variables under CONFIGURATION.

Tested and working using:
'EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml'

EOF
}

# ----------------- Configuration ----------------------------------------

# Species specific sample directory, abosolute path
SD=

# Working directory, abosolute path
WD="$(readlink -f "$PWD")"

# Algorithm to use during alignment: mem (>70MB, contemporary samples) [default] or aln (<70MB, historic samples)
algo="mem"

# Define memory per cpu in G (must be integer)
memory="8"

# Define number of cpus to be used (must be integer)
cpus="8"

# Option to force run even if target sample directory contains .bam
force_overwrite="N"

# Name of data preparation directory
dataprep="03_init_analysis_files"

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

# ----------------- Script Flag Processing -------------------------------

while getopts 's:d:m:c:fh' OPTION; do
    case "$OPTION" in
        s)
            if [ -d "$OPTARG" ]; then
                SD="$(readlink -f "$OPTARG")"
            else
                echo -e "\nERROR: $OPTARG does not seem to be a directory\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
                exit 1
            fi
            ;;
        d)
            WD="$(readlink -f "$OPTARG")"
            ;;
        m)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                memory="$OPTARG"
            else
                echo -e "\nERROR: -m must be a whole integer, 8 [default]'\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
                exit 1
            fi
            ;;
        c)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                cpus="$OPTARG"
            else
                echo -e "\nERROR: -c must be a whole integer, 8 [default]'\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
                exit 1
            fi
            ;;
        f)
            force_overwrite="Y"
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

# Removes "/" from the end of path for species directory
if [ "${SD: -1}" == "/" ]; then
    length=${#SD}
    SD=${SD::length - 1}
fi

# Removes "/" from the end of path for working directory
if [ "${WD: -1}" == "/" ]; then
    length=${#WD}
    WD=${WD::length - 1}
fi

# Tests that species directory has been supplied
if [[ -z $SD ]]; then
    usage
    exit 1
fi

# ----------------- Functions --------------------------------------------

# Function used to adjust time for jobs
timer()
{
    awk -v adjustment="$1" -v base="$2" -v static="$3" 'BEGIN { print int( base * adjustment + static) }'
}

# Function used to queue jobs
queue() # "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo" "$memory" "$adjustment" "$script_path" "$age"
{
    # Location of logfile
    logfile="$4"/"$6"/"$(basename "$3")"/log_"$(basename "$3")"_"$(date +"%Y-%m-%d")".txt

    stdoutput="$4"/"$6"/"$(basename "$3")"/"$(basename "$5")"/stdout_sbatch

    echo "Initial Analysis $(basename "$5") is sent to the queue" >> "$logfile"

    # Mpileup
    jid1=$(sbatch \
            --parsable \
            --time="$(timer "$9" 720 120)" \ # ***Figure out time***
            --mem-per-cpu="$8"G \
            --cpus-per-task="$1" \
            --output="$stdoutput"/"$(basename "$5")"_01-%j.out \
            "${10}"/modules/03_01_mpileup.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7") # ***Add actual arguments***
                            
    echo -e "\t'Mpileup' job has been submitted for $(basename "$5") -- Job ID: $jid1" >> "$logfile"

    # VCF
    jid2=$(sbatch \
            --parsable \
            --time="$(timer "$9" 1800 120)" \ # ***Figure out time***
            --mem-per-cpu="$8"G \
            --cpus-per-task="$1" \
            --output="$stdoutput"/"$(basename "$5")"_02-%j.out \
            --dependency=afterany:"$jid1" \
            "${10}"/modules/03_02_vcf.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7") # ***Add actual arguments***

    echo -e "\t'VCF' job has been submitted for $(basename "$5") -- Job ID: $jid2" >> "$logfile"

    # Split file
    jid3=$(sbatch \
            --parsable \
            --time="$(timer "$9" 1800 120)" \ # ***Figure out time***
            --mem-per-cpu="$8"G \
            --cpus-per-task="$1" \
            --output="$stdoutput"/"$(basename "$5")"_03-%j.out \
            --dependency=afterany:"$jid1" \
            "${10}"/modules/03_03_split.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7") # ***Add actual arguments***

    echo -e "\t'File Split' job has been submitted for $(basename "$5") -- Job ID: $jid3" >> "$logfile"

    # SFS queue looper
    jid4=$(sbatch \
            --parsable \
            --time="$(timer "$9" 1800 120)" \ # ***Figure out time***
            --mem-per-cpu="$8"G \
            --cpus-per-task="$1" \
            --output="$stdoutput"/"$(basename "$5")"_04-%j.out \
            --dependency=afterany:"$jid3" \
            "${10}"/modules/03_04_sfs_looper.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7") # ***Add actual script name*** ***Add actual arguments***

    echo -e "\t'SFS Looper' job has been submitted for $(basename "$5") -- Job ID: $jid4" >> "$logfile"

    mapfile -t id < <(squeue -u "$USER" -S -V -o %A)

    echo "${id[1]}"

    # Clean up of empty stdout files
    sbatch \
            --output=/dev/null \
            --error=/dev/null \
            --dependency=afterany:"${id[1]}" \
            "${10}"/modules/02_08_cleanup.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7"
}

# ----------------- Script Queue -----------------------------------------

# Holds path to script directory
script_path=$(dirname "$script_path")

# Finds references genome for designated species or exits if it cannot be found
for ref in /faststorage/project/EcoGenetics/BACKUP/reference_genomes/"$(basename "$SD")"/*.fna; do
    if [ ! -e "$ref" ]; then
        echo -e "\nCannot locate reference genome, $ref\n"
        exit 1
    else
        RG=$ref
        break
    fi
done

# Exit with error message if reference genome is not indexed
for index in "$(dirname "$RG")"/*.ann; do
    if [ ! -e "$index" ]; then
        echo -e "\nDesignated reference genome does not seem to be indexed, $RG\n"
        exit 1
    fi
done

# Creates working directory if it doesn't exist
[[ -d "$WD" ]] || mkdir -m 775 "$WD"

# Creates temp directory in working directory if none exist
[[ -d "$WD"/temp ]] || mkdir -m 775 "$WD"/temp

# Creates data preparation directory in working directory if none exist
[[ -d "$WD"/"$dataprep" ]] || mkdir -m 775 "$WD"/"$dataprep"

# Creates species directory in data_preparation directory if none exist
speciesdir="$WD"/"$dataprep"/"$(basename "$SD")"
[[ -d "$speciesdir" ]] || mkdir -m 775 "$speciesdir"

# Creates log file
touch "$speciesdir"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
logfile="$speciesdir"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

# Makes sure all associated modules are executable
for script in "$script_path"/modules/*; do
    [ ! -x "$script" ] && chmod -R u+rwx "$script_path"/modules && break
done

# Loops through all sample folders within species specific sample directory
for sample in "$SD"/*; do

        # Checks if sample folder is empty
        if [ "$(ls -A "$sample")" ]; then
        
                # Checks whether a .bam file already exists within sample folder, indicating samples have already been processed        
                for file in "$sample"/*.bam; do
                        if [ ! -e "$file" ] || [ "$force_overwrite" == "Y" ]; then
                
                                # Checks sample file size
                                filesize=
                                for fna in "$sample"/*.fq.gz; do
                                        [ ! "$filesize" ] && filesize=$(wc -c < "$fna") && break
                                done

                                # Change requested time on nodes depending on filesize comparative to R1 Ocin_NYS-F
                                adjustment=$(awk -v filesize="$filesize" 'BEGIN { print ( filesize / 93635424798 + 0.1) }')
                                
                                sampledir="$speciesdir"/"$(basename "$sample")"

                                # Creates sample directory in species directory if none exist
                                [[ -d "$sampledir" ]] || mkdir -m 775 "$sampledir"

                                # Creates pre- and post-filtering directory in sample directory if none exist
                                [[ -d "$sampledir"/pre_filter_stats ]] || mkdir -m 775 "$sampledir"/pre_filter_stats
                                [[ -d "$sampledir"/post_filter_stats ]] || mkdir -m 775 "$sampledir"/post_filter_stats

                                # Create folder for STDOUT files generated sbatch
                                [[ -d "$sampledir"/stdout_sbatch ]] || mkdir -m 775 "$sampledir"/stdout_sbatch

                                        # Checks if currently working with museomics samples
                                        if [[ "$SD" == *"museomics"* ]]; then

                                                # Checks if sample directory is pre-2000 (historic) or post-2000 (modern) and sets the algrorithm parameter correspondingly
                                                if [ "$((${sample: -4}))" -gt 2000 ]; then

                                                        # Sample is modern
                                                        age=" modern"
                                                        algo="mem"
                        
                                                        queue "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo" "$memory" "$adjustment" "$script_path" "$age"

                                                else

                                                        # Sample is historic
                                                        age=" historic"
                                                        algo="aln"

                                                        queue "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo" "$memory" "$adjustment" "$script_path" "$age"
                        
                                                fi

                                        else

                                                # Sample is custom
                                                age=""

                                                queue "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo" "$memory" "$adjustment" "$script_path" "$age"
                                        fi

                        else

                                echo -e "$(basename "$sample") already contains a .bam file, $(basename "$file"), and is skipped\n" >> "$logfile"

                        fi

                done

        else

                echo -e "$(basename "$sample") is an empty directory and is skipped\n" >> "$logfile"

        fi

done

echo -e "\nLog file can be found here: $logfile\n"

exit 0
