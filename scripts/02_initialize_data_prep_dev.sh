#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:40:00
#SBATCH --output=Data_Prep_Initialize-%j.out

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

Usage: 02_initialize_data_prep.sh [PARAMETERS] [OPTIONS]

This script is used for initializing the standardized data preparation procedure for sequence data.
Intended to be used in conjunction with 'sbatch', however it also be used in conjunction with 'srun' 
or simply on the frontend as its resource-demands are fairly low.

The 'Mads loop condition': If the specified species sample directory is within the 'museomics' directory, 
the script will automatically try to detect whether the sample directory is contemporary or historical 
and automatically applies the corresponding algorithm. This means that if your samples are in the 
museomics directory you do not need to specify an algorithm.

PARAMETERS (must be assigned):
    -r  FILE            Species specific reference genome, abosolute path (reference genome in FASTA format)
    -s  DIRECTORY       Species specific sample directory, abosolute path (Do NOT end with '/')
    -d  DIRECTORY       Working directory, abosolute path (Do NOT end with '/')

OPTIONS:
    -a  ALGORITHM       Choice of algorithm to be used during alignment. 'mem' (>70MB, contemporary samples)[default] or 'aln' (<70MB, historic samples)
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
RG=

# Species specific sample directory, abosolute path (Do NOT end with '/')
SD=

# Working directory, abosolute path (Do NOT end with '/')
WD=

# Algorithm to use during alignment: mem (>70MB, contemporary samples) [default] or aln (<70MB, historic samples)
algo="mem"

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

error_a()
{
cat << EOF

ERROR: -a must be either 'mem' [default] or 'aln'

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

while getopts 'r:s:d:a:m:c:h' OPTION; do
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
[[ -d "$WD" ]] || mkdir -m 775 "$WD" # || mkdir -m 764 "$PWD"/"$WD"

# Creates temp directory in working directory if none exist
[[ -d "$WD"/temp ]] || mkdir -m 775 "$WD"/temp

# Creates data_preparation directory in working directory if none exist
[[ -d "$WD"/"$(basename "$script_path")" ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"

# Creates species directory in data_preparation directory if none exist
[[ -d "$WD"/"$(basename "$script_path")"/"$(basename "$SD")" ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"

# Creates log file
touch "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

# Loops through all sample folders within species specific sample directory
for sample in "$SD"/*; do

    variables=("$RG" "$SD" "$WD" "$algo" "$memory" "$cpus" "$sample")

    # Checks if sample folder is empty
    if [ "$(ls -A "$sample")" ]; then
        
        # Checks whether a .bam file already exists within sample folder, indicating samples have already been processed        
        for file in "$sample"/*.bam; do
            if [ ! -e "$file" ]; then
                
                # Checks sample file size
                filesize=
                for fna in "$sample"/*.fna; do
                    if [ ! "$filesize" ]; then
                        filesize=$(wc -c < "$fna")
                    fi
                done

                # Change requested time on nodes depending on filesize comparative to R1 Ocin_NYS-F
                adjustment=$(awk -v filesize="$filesize" 'BEGIN { print ( filesize / 93635424798 ) }')

                # Creates sample directory in species directory if none exist
                [[ -d "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")" ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"

                # Creates pre- and post-filtering directory in sample directory if none exist
                [[ -d "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats
                [[ -d "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats

                # Create folder for STDOUT files generated sbatch
                [[ -d "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch ]] || mkdir -m 775 "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch

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
                            
                            echo "$(basename "$sample") is sent to the queue as a single-end modern sample" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                            # AdapterRemoval
                            jid1=$(sbatch \
                                    --parsable \
                                    --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 720 * adjustment + 120) }')" \
                                    --mem-per-cpu="$memory"G \
                                    --cpus-per-task="$cpus" \
                                    --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_01-%j.out \
                                    "$script_path"/modules/02_01_single_adapterremoval.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                            
                            echo "'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: $jid1" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                            # Aligning to reference
                            jid2=$(sbatch \
                                    --parsable \
                                    --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 1800 * adjustment + 120) }')" \
                                    --mem-per-cpu="$memory"G \
                                    --cpus-per-task="$cpus" \
                                    --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_02-%j.out \
                                    --dependency=afterany:"$jid1" \
                                    "$script_path"/modules/02_02_single_alignment.sh "$RG" "$SD" "$WD" "$sample" "$algo" "$cpus")

                            echo "'Alignment' job has been submitted for $(basename "$sample") -- Job ID: $jid2" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                            # Marking duplicates
                            jid4=$(sbatch \
                                    --parsable \
                                    --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 1800 * adjustment + 120) }')" \
                                    --mem-per-cpu="$memory"G \
                                    --cpus-per-task="$cpus" \
                                    --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_04-%j.out \
                                    --dependency=afterany:"$jid2" \
                                    "$script_path"/modules/02_04_markduplicates.sh "$RG" "$SD" "$WD" "$sample" "$cpus")

                            echo "'Marking duplicates' job has been submitted for $(basename "$sample") -- Job ID: $jid4" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        
                        else

                            # Paired-end
                            
                            echo "$(basename "$sample") is sent to the queue as a pair-end modern sample" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                            # AdapterRemoval
                            jid1=$(sbatch \
                                    --parsable \
                                    --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 720 * adjustment +120 ) }')" \
                                    --mem-per-cpu="$memory"G \
                                    --cpus-per-task="$cpus" \
                                    --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_01-%j.out \
                                    "$script_path"/modules/02_01_paired_adapterremoval.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                            
                            echo "'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: $jid1" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                            # Aligning to reference
                            jid2_1=$(sbatch \
                                    --parsable \
                                    --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 1800 * adjustment +120 ) }')" \
                                    --mem-per-cpu="$memory"G \
                                    --cpus-per-task="$cpus" \
                                    --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_02_01-%j.out \
                                    --dependency=afterany:"$jid1" \
                                    "$script_path"/modules/02_02_paired_01_alignment.sh "$RG" "$SD" "$WD" "$sample" "$algo" "$cpus")
                            jid2_2=$(sbatch \
                                    --parsable \
                                    --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 1800 * adjustment +120 ) }')" \
                                    --mem-per-cpu="$memory"G \
                                    --cpus-per-task="$cpus" \
                                    --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_02_02-%j.out \
                                    --dependency=afterany:"$jid1" \
                                    "$script_path"/modules/02_02_paired_02_alignment.sh "$RG" "$SD" "$WD" "$sample" "$algo" "$cpus")

                            echo "'Alignment of paired ends' job has been submitted for $(basename "$sample") -- Job ID: $jid2_1" >> "$WD"/0"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                            echo "'Alignment of collapsed single end' job has been submitted for $(basename "$sample") -- Job ID: $jid2_2" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                            # Merging of alignment files
                            jid3=$(sbatch \
                                    --parsable \
                                    --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 1200 * adjustment +120 ) }')" \
                                    --mem-per-cpu="$memory"G \
                                    --cpus-per-task="$cpus" \
                                    --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_03-%j.out \
                                    --dependency=afterany:"$jid2_1":"$jid2_2" \
                                    "$script_path"/modules/02_03_paired_merge.sh "$RG" "$SD" "$WD" "$sample" "$cpus")

                            echo "'Merging' job has been submitted for $(basename "$sample") -- Job ID: $jid3" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                            # Marking duplicates
                            jid4=$(sbatch \
                                    --parsable \
                                    --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 1800 * adjustment +120 ) }')" \
                                    --mem-per-cpu="$memory"G \
                                    --cpus-per-task="$cpus" \
                                    --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_04-%j.out \
                                    --dependency=afterany:"$jid3" \
                                    "$script_path"/modules/02_04_markduplicates.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                            
                            echo "'Marking duplicates' job has been submitted for $(basename "$sample") -- Job ID: $jid4" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                        fi
                        
                        # Statistics pre-filtering
                        jid5_1=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 60 * adjustment + 60) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_05_01-%j.out \
                                --dependency=afterany:"$jid4" \
                                "$script_path"/modules/02_05_01_prestats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        jid5_2=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 60 * adjustment + 60) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_05_02-%j.out \
                                --dependency=afterany:"$jid4" \
                                "$script_path"/modules/02_05_02_prestats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        jid5_3=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 360 * adjustment + 120 ) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_05_03-%j.out \
                                --dependency=afterany:"$jid4" \
                                "$script_path"/modules/02_05_03_prestats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        jid5_4=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 60 * adjustment + 60) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_05_04-%j.out \
                                --dependency=afterany:"$jid4" \
                                "$script_path"/modules/02_05_04_prestats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        
                        echo "'Pre-filtering statistics (idxstats)' job has been submitted for $(basename "$sample") -- Job ID: $jid5_1" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        echo "'Pre-filtering statistics (flagstat)' job has been submitted for $(basename "$sample") -- Job ID: $jid5_2" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        echo "'Pre-filtering statistics (coverage)' job has been submitted for $(basename "$sample") -- Job ID: $jid5_3" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        echo "'Pre-filtering statistics (stats)' job has been submitted for $(basename "$sample") -- Job ID: $jid5_4" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt


                        # Removal of duplicates, unmapped reads and low quality mappings
                        jid6=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 1440 * adjustment + 120 ) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_06-%j.out \
                                --dependency=afterany:"$jid5_1":"$jid5_2":"$jid5_3":"$jid5_4" \
                                "$script_path"/modules/02_06_filtering.sh "$RG" "$SD" "$WD" "$sample" "$cpus")

                        echo "'Filtering' job has been submitted for $(basename "$sample") -- Job ID: $jid6" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                        # Statistics post-filtering
                        jid7_1=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 60 * adjustment + 60 ) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_07_01-%j.out \
                                --dependency=afterany:"$jid6" \
                                "$script_path"/modules/02_07_01_poststats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        jid7_2=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 60 * adjustment + 60 ) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_07_02-%j.out \
                                --dependency=afterany:"$jid6" \
                                "$script_path"/modules/02_07_02_poststats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        jid7_3=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 360 * adjustment + 120 ) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_07_03-%j.out \
                                --dependency=afterany:"$jid6" \
                                "$script_path"/modules/02_07_03_poststats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        jid7_4=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 360 * adjustment + 120 ) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_07_04-%j.out \
                                --dependency=afterany:"$jid6" \
                                "$script_path"/modules/02_07_04_poststats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        jid7_5=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 60 * adjustment + 60 ) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_07_05-%j.out \
                                --dependency=afterany:"$jid6" \
                                "$script_path"/modules/02_07_05_poststats.sh "$RG" "$SD" "$WD" "$sample" "$cpus")
                        jid7_6=$(sbatch \
                                --parsable \
                                --time="$(awk -v adjustment="$adjustment" 'BEGIN { print int( 1200 * adjustment + 120 ) }')" \
                                --mem-per-cpu="$memory"G \
                                --cpus-per-task="$cpus" \
                                --output="$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch/"$(basename "$sample")"_07_06-%j.out \
                                --dependency=afterany:"$jid6" \
                                "$script_path"/modules/02_07_06_qualimap.sh "$RG" "$SD" "$WD" "$sample" "$cpus")

                        echo "'Post-filtering statistics (flagstat)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_1" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        echo "'Post-filtering statistics (idxstats)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_2" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        echo "'Post-filtering statistics (coverage)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_3" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        echo "'Post-filtering statistics (readchange)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_4" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        echo "'Post-filtering statistics (stats)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_5" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt
                        echo "'Post-filtering statistics (qualimap)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_6" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                        # Clean up of empty stdout files
                        sbatch \
                                --output=/dev/null \
                                --error=/dev/null \
                                --dependency=afterany:"$jid7_1":"$jid7_2":"$jid7_3":"$jid7_4":"$jid7_5":"$jid7_6" \
                                "$script_path"/modules/02_08_data_prep.sh "$RG" "$SD" "$WD" "$sample"

                    else

                        # Sample is historic
                        algo="aln"

                        # Check the number of .fq.gz files in sample directory (assumed to be indicative of whether sample is single- or paired-end)
                        count=$(find "$sample"/ -maxdepth 1 -type f -name '*.fq.gz' | wc -l)
                        if [ "$count" == 1 ]; then
                            
                            # Single-end

                            echo "$(basename "$sample") is sent to the queue as a single-end historic sample" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                        else

                            # Paired-end
                            
                            echo "$(basename "$sample") is sent to the queue as a pair-end historic sample" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                        fi

                    fi

                else

                    # Check the number of .fq.gz files in sample directory (assumed to be indicative of whether sample is single- or paired-end)
                    count=$(find "$sample"/ -maxdepth 1 -type f -name '*.fq.gz' | wc -l)
                    if [ "$count" == 1 ]; then

                        # Single-end
                        
                        echo "$(basename "$sample") is sent to the queue as a single-end sample" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                    else

                        # Paired-end

                        echo "$(basename "$sample") is sent to the queue as a pair-end sample" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

                    fi

                fi

            else

                echo "$(basename "$sample") already contains a .bam file, $(basename "$file"), and is skipped" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

            fi

        done

    else

        echo "$(basename "$sample") is an empty directory and is skipped" >> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

    fi

done

echo "Log file can be found here: $WD/$(basename "$script_path")/$(basename "$SD")/log_$(basename "$SD")_$(date +"%Y-%m-%d").txt"

exit 0

