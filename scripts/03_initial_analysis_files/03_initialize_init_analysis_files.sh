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

# ----------------- Environment ------------------------------------------

if [ "$USER" == "jepe" ]; then
    source /home/"$USER"/.bashrc
    source activate ecogen_primary
fi

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
    -a  TYPE            Choice of analysis files to create. 'full' = the entire genome,
                        'intergene' = for intergenic regions, 'nonsyn' = for non-
                        synonymous regions, 'all' for all of the above [default]
    -m  INTEGER         Amount of memory to be used by each CPU. 8 [default]
    -c  INTEGER         Number of CPUs to be used. 8 [default]
    -u                  Run only on a single sample. -s then needs to lead to
						the specific sample
    -f                  Force run even if target sample directory contains relevant files
    -h                  Show this message

Parameters can also be set directly in the script. 
If you decide to do so, change the variables under CONFIGURATION.

Tested and working using:
'EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml'

EOF
}

# ----------------- Configuration ----------------------------------------

# Species specific sample directory
SD=

# Working directory, abosolute path
WD="$(readlink -f "$PWD")"

# Analysis files to create
queuetype="all"

# Number of parts to split pileup into
parts="100"

# Define memory per cpu in G (must be integer)
memory="8"

# Define number of cpus to be used (must be integer)
cpus="8"

# Option to only run a single sample
single_sample="N"

# Option to force run even if target sample directory contains .bam
force_overwrite="N"

# Name of data directory
data="03_init_analysis_files"

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

# ----------------- Script Flag Processing -------------------------------

while getopts 's:d:a:p:m:c:ufh' OPTION; do
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
        a)
            if [ "$OPTARG" == "all" ] || [ "$OPTARG" == "full" ]  || [ "$OPTARG" == "intergene" ] || [ "$OPTARG" == "nonsyn" ]; then
                queuetype="$OPTARG"
            else
                echo -e "\nERROR: $OPTARG is not a recognized choice\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
                exit 1
            fi
            ;;
        p)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                parts="$OPTARG"
            else
                echo -e "\nERROR: -p must be a whole integer, 20 [default]'\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
                exit 1
            fi
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
        u)
			single_sample="Y"
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

# Single sample argument correction
if [ $single_sample == "Y" ]; then
	sample="$SD"
	SD=$(dirname "$SD")
fi

# ----------------- Functions --------------------------------------------

# Function for sample processing
sample_processing()
{
	# Checks if sample folder is empty
	if [ "$(ls -A "$sample")" ]; then
		
		# Checks if directory contains .bam file       
		for bam in "$sample"/*.bam; do
			if [ -e "$bam" ]; then

                # # Checks sample file size
                # filesize=
                # for fna in "$sample"/*.fq.gz; do
                # 	[ ! "$filesize" ] && filesize=$(wc -c < "$fna") && break
                # done

                # # Change requested time on nodes depending on filesize comparative to R1 Ocin_NYS-F
                # adjustment=$(awk -v filesize="$filesize" 'BEGIN { print ( filesize / 93635424798 + 0.1) }')

                # Creates sample directory in species directory if none exist
                sampledir="$speciesdir"/"$(basename "$sample")"
                [[ -d "$sampledir" ]] || mkdir -m 775 "$sampledir"

                # Create folder for STDOUT files generated sbatch
                stdoutput="$sampledir"/stdout_sbatch
                [[ -d "$stdoutput" ]] || mkdir -m 775 "$stdoutput"

                # Queues jobs for samples
                if [ "$queuetype" == "all" ]; then

                    echo "Initial Analysis Files (All) for $(basename "$sample") is sent to the queue" >> "$logfile"

                    queue01

                    queue02

                    queue03
                
                elif [ "$queuetype" == "full" ]; then

                    echo "Initial Analysis Files (Full) for $(basename "$sample") is sent to the queue" >> "$logfile"

                    queue01

                elif [ "$queuetype" == "intergene" ]; then

                    echo "Initial Analysis Files (Intergenic) for $(basename "$sample") is sent to the queue" >> "$logfile"

                    queue02
                
                elif [ "$queuetype" == "nonsyn" ]; then

                    echo "Initial Analysis Files (Non-synonymous) for $(basename "$sample") is sent to the queue" >> "$logfile"

                    queue03
                
                fi

			else

				echo -e "$(basename "$sample") doesn't contain a .bam file, $(basename "$bam"), and is skipped\n" >> "$logfile"

			fi

		done

	else

		echo -e "$(basename "$sample") is an empty directory and is skipped\n" >> "$logfile"

	fi
}

# Function used to adjust time for jobs
timer()
{
    awk -v adjustment="$1" -v base="$2" -v static="$3" 'BEGIN { print int( base * adjustment + static) }'
}

# Functions used to queue jobs
queue01()
{
    if [ ! -e "$sample"/"$(basename "$sample")".pileup ]; then

        # Create pileup
        jid1_1=$(sbatch \
            --parsable \
            --time=1440 \
            --mem-per-cpu="$memory"G \
            --cpus-per-task="$cpus" \
            --output="$stdoutput"/"$(basename "$sample")"_01_01-%j.out \
            "$script_path"/modules/03_01_01_full_pileup.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data") # ***Add actual arguments***
                                
        echo -e "\t'Pileup - Full' job has been submitted for $(basename "$sample") -- Job ID: $jid1_1" >> "$logfile"

        # Split pileup
        jid2_1=$(sbatch \
            --parsable \
            --time=1440 \
            --mem-per-cpu="$memory"G \
            --cpus-per-task="$cpus" \
            --output="$stdoutput"/"$(basename "$sample")"_02_01-%j.out \
            --dependency=afterany:"$jid1_1" \
            "$script_path"/modules/03_02_01_split_full.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data") # ***Add actual arguments***

        echo -e "\t'Split Pileup - Full' job has been submitted for $(basename "$sample") -- Job ID: $jid2_1" >> "$logfile"
    
    else

        # Split pileup
        jid2_1=$(sbatch \
            --parsable \
            --time=1440 \
            --mem-per-cpu="$memory"G \
            --cpus-per-task="$cpus" \
            --output="$stdoutput"/"$(basename "$sample")"_02_01-%j.out \
            "$script_path"/modules/03_02_01_split_full.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data") # ***Add actual arguments***
        
        echo -e "\tUsing existing full pileup file for $(basename "$sample")" >> "$logfile"
        echo -e "\t'Split Pileup - Full' job has been submitted for $(basename "$sample") -- Job ID: $jid2_1" >> "$logfile"

    fi

    # Full SFS array
    jid3_1=$(sbatch \
        --parsable \
        --array=0-"$parts" \
        --chdir="$WD"/temp \
        --time=1440 \
        --mem-per-cpu="$memory"G \
        --cpus-per-task="$cpus" \
        --output="$stdoutput"/"$(basename "$sample")"_03_01-%j.out \
        --dependency=afterany:"$jid2_1" \
        "$script_path"/modules/03_03_01_sfs_full.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data") # ***Add actual arguments***

    echo -e "\t'Site Ferquency Spectrum - Full' job has been submitted for $(basename "$sample") -- Job ID: $jid3_1" >> "$logfile"

    # Full SFS assembly
    # jid4_1=$(sbatch \
    #   --parsable \
    #   --time=1440 \
    #   --mem-per-cpu="$memory"G \
    #   --cpus-per-task="$cpus" \
    #   --output="$stdoutput"/"$(basename "$sample")"_04_01-%j.out \
    #   --dependency=afterany:"$jid3_1" \
    #   "$script_path"/modules/03_04_01_full_assemble.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data") # ***Add actual arguments***

    # echo -e "\t'Site Frequency Spectrum Assembly - Full' job has been submitted for $(basename "$sample") -- Job ID: $jid4_1\n" >> "$logfile"

    # Clean up of empty stdout files
    # sbatch \
    #   --output=/dev/null \
    #   --error=/dev/null \
    #   --dependency=afterany:"${id[1]}" \
    #   "$script_path"/modules/02_08_cleanup.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data"
}

queue02()
{
    if [ ! -e "$sampledir"/"$(basename "$sample")"_intergenic.pileup ]; then
        
        # Create pileup
        jid1_2=$(sbatch \
            --parsable \
            --time=1440 \
            --mem-per-cpu="$memory"G \
            --cpus-per-task="$cpus" \
            --output="$stdoutput"/"$(basename "$sample")"_01_02-%j.out \
            "$script_path"/modules/03_01_02_intergenic_pileup.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***
                                
        echo -e "\t'Pileup - Intergenic' job has been submitted for $(basename "$sample") -- Job ID: $jid1_2" >> "$logfile"

        # Split pileup
        jid2_2=$(sbatch \
            --parsable \
            --time=1440 \
            --mem-per-cpu="$memory"G \
            --cpus-per-task="$cpus" \
            --output="$stdoutput"/"$(basename "$sample")"_02_02-%j.out \
            --dependency=afterany:"$jid1_2" \
            "$script_path"/modules/03_02_02_split_intergenic.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***

        echo -e "\t'Split Pileup - Intergenic' job has been submitted for $(basename "$sample") -- Job ID: $jid2_2" >> "$logfile"

    else

        # Split pileup
        jid2_2=$(sbatch \
            --parsable \
            --time=1440 \
            --mem-per-cpu="$memory"G \
            --cpus-per-task="$cpus" \
            --output="$stdoutput"/"$(basename "$sample")"_02_02-%j.out \
            "$script_path"/modules/03_02_02_split_intergenic.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***

        echo -e "\tUsing existing intergenic pileup file for $(basename "$sample")" >> "$logfile"
        echo -e "\t'Split Pileup - Intergenic' job has been submitted for $(basename "$sample") -- Job ID: $jid2_2" >> "$logfile"
    
    fi

    # SFS array
    jid3_2=$(sbatch \
        --parsable \
        --array=0-"$parts" \
        --chdir="$WD"/temp \
        --time=1440 \
        --mem-per-cpu="$memory"G \
        --cpus-per-task="$cpus" \
        --output="$stdoutput"/"$(basename "$sample")"_03_02-%j.out \
        --dependency=afterany:"$jid2_2" \
        "$script_path"/modules/03_03_02_sfs_intergenic.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***

    echo -e "\t'Site Ferquency Spectrum - Intergenic' job has been submitted for $(basename "$sample") -- Job ID: $jid3_2" >> "$logfile"

    # SFS assembly
    # jid4_2=$(sbatch \
    #   --parsable \
    #   --time=1440 \
    #   --mem-per-cpu="$memory"G \
    #   --cpus-per-task="$cpus" \
    #   --output="$stdoutput"/"$(basename "$sample")"_04_02-%j.out \
    #   --dependency=afterany:"$jid3_2" \
    #   "$script_path"/modules/03_04_02_intergenic_assemble.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***

    # echo -e "\t'Site Frequency Spectrum Assembly - Intergenic' job has been submitted for $(basename "$sample") -- Job ID: $jid4_2\n" >> "$logfile"

    # Clean up of empty stdout files
    # sbatch \
    #   --output=/dev/null \
    #   --error=/dev/null \
    #   --dependency=afterany:"${id[1]}" \
    #   "$script_path"/modules/02_08_cleanup.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data "$parts" "$script_path""
}

queue03()
{
    # Create pileup
    jid1_3=$(sbatch \
        --parsable \
        --time=1440 \
        --mem-per-cpu="$memory"G \
        --cpus-per-task="$cpus" \
        --output="$stdoutput"/"$(basename "$sample")"_01_03-%j.out \
        "$script_path"/modules/03_01_03*****.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***
                            
    echo -e "\t'Pileup - Non-synonymous' job has been submitted for $(basename "$sample") -- Job ID: $jid1_3" >> "$logfile"

    # Split pileup
    jid2_3=$(sbatch \
        --parsable \
        --time=1440 \
        --mem-per-cpu="$memory"G \
        --cpus-per-task="$cpus" \
        --output="$stdoutput"/"$(basename "$sample")"_02_03-%j.out \
        --dependency=afterany:"$jid1_3" \
        "$script_path"/modules/03_02_03*****.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***

    echo -e "\t'Split Pileup - Non-synonymous' job has been submitted for $(basename "$sample") -- Job ID: $jid2_3" >> "$logfile"

    #SFS array
    jid3_3=$(sbatch \
        --parsable \
        --array=0-"$parts" \
        --chdir="$WD"/temp \
        --time=1440 \
        --mem-per-cpu="$memory"G \
        --cpus-per-task="$cpus" \
        --output="$stdoutput"/"$(basename "$sample")"_03_03-%j.out \
        --dependency=afterany:"$jid2_3" \
        "$script_path"/modules/03_03_03*****.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***

    echo -e "\t'Site Ferquency Spectrum - Non-synonymous' job has been submitted for $(basename "$sample") -- Job ID: $jid3_3" >> "$logfile"

    # Full SFS assembly
    jid4_3=$(sbatch \
        --parsable \
        --time=1440 \
        --mem-per-cpu="$memory"G \
        --cpus-per-task="$cpus" \
        --output="$stdoutput"/"$(basename "$sample")"_04_03-%j.out \
        --dependency=afterany:"$jid3_3" \
        "$script_path"/modules/03_04_03*****.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data" "$parts" "$script_path") # ***Add actual arguments***

    echo -e "\t'Site Frequency Spectrum Assembly - Non-synonymous' job has been submitted for $(basename "$sample") -- Job ID: $jid4_3" >> "$logfile"

    # Clean up of empty stdout files
    # sbatch \
    #   --output=/dev/null \
    #   --error=/dev/null \
    #   --dependency=afterany:"${id[1]}" \
    #   "$script_path"/modules/02_08_cleanup.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$data "$parts" "$script_path""
}

# ----------------- Script Queue -----------------------------------------

# Adjusts number of parts as count starts with 0
parts=$(awk -v p="$parts" 'BEGIN { print ( p - 1 ) }')

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

# Creates data directory in working directory if none exist
[[ -d "$WD"/"$data" ]] || mkdir -m 775 "$WD"/"$data"

# Creates species directory in data directory if none exist
speciesdir="$WD"/"$data"/"$(basename "$SD")"
[[ -d "$speciesdir" ]] || mkdir -m 775 "$speciesdir"

# Creates log file
touch "$speciesdir"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
logfile="$speciesdir"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

# Makes sure all associated modules are executable
for script in "$script_path"/modules/*; do
    [ ! -x "$script" ] && chmod -R u+rwx "$script_path"/modules && break
done

# Checks whether or not to run on single sample
if [ "$single_sample" == "N" ]; then

	# Multiple samples

	# Loops through all sample folders within species specific sample directory
	for sample in "$SD"/*; do

		sample_processing
		
	done

else

	# Single sample
	sample_processing

fi

echo -e "\nLog file can be found here: $logfile\n"

exit 0
