#!/bin/bash

# ------------------ Description ------------------------------------------

# Script for initializing data preparation procedure for re-sequencing data

# ----------------- Usage and version -------------------------------------

usage(){
cat << EOF

Usage: 02_initialize_data_prep.sh -s <FILE> [OPTIONS]

This script is used for initializing the standardized data preparation procedure
for sequence data at the Center for EcoGenetics.

If the specified species sample directory is within the 'museomics' directory,
the script will automatically try to detect whether the sample directory is
contemporary or historical and automatically apply the corresponding algorithm.
This means that if your samples are in the museomics directory you do not need
to specify an algorithm.

    -s  DIRECTORY           Species specific directory.

OPTIONS:
    -d  DIRECTORY           Working directory. If not assigned, will use current 
                            working directory [default].
    -a  ALGORITHM           Choice of algorithm to be used during alignment.
                            'mem' (>70MB, contemporary samples)[default] or 
                            'aln' (<70MB, historic samples).
    -r  FILE                Manually set reference genome.
    -u                      Run only on a single sample. -s then needs to lead to
                            the specific sample.
    -f                      Force run even if target sample directory contains .bam.
    -j  STRING              If 'list' shows jobgraph. Can be supplied with a single
                            integer or a range of integer to define which jobs in the
                            pipeline should be run.
    -t  INTEGER             Timefactor.
    -h                      Show this message.
    -v                      Show current version and rescent version history.

EOF
}

# Current version and rescent version history
version()
{
cat << EOF

Current Version: 1.3.1
Updated 18/04-2023

Version history:
    Version             Changelog
    1.0                 Added version information and updated
                        timer for samtools stats function and script 05_01.
	1.1					If sample is marked as 'historic' more lenient
						criteria will be used during alignment.
	1.1.1				Temporary fix to usage of 'aln/samse/sampe'
						algorithm.
	1.2.0				Added script for getting read depth post filtering
						and script for extracted unmapped read. Other minor
						changes have been made as well.
	1.2.1				Removes .out and log files left over from previous
						runs.
    1.3.0               Major structural overhaul to code.
    1.3.1               Added posibility to manually choose reference genome.

EOF
}

# ----------------- Configuration ----------------------------------------

# Species specific sample directory, abosolute path
SD=

# Working directory
WD="$(readlink -f "$PWD")"

# Algorithm to use during alignment: mem (>70MB, contemporary samples) [default] or aln (<70MB, historic samples)
algo="mem"

# Option to only run a single sample
single_sample="N"

# Option to force run even if target sample directory contains .bam
force_overwrite="N"

# Name of created data directory under the working directory
datadir="02_data_preparation"

# Timefactor
timefactor=1

# Makes extended globbing available
shopt -s extglob

# Creates list for the execution of jobs
num_of_jobs=23
for (( i=1; i<="$num_of_jobs"; i++ )); do
    joblist["$i"]=true
done

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command={print $2}')
# If run on the frontend:
else
    script=$(realpath "$0")
fi
script_path="$(dirname "$script")"

jobgraph(){
cat << EOF

1. AdapterRemoval (single-end)      2. AdapterRemoval (paired-end)      3. AdapterRemoval (several paired-end)
            | |                                  \ \___________________________________/ /
            | |                                   \_________________   _________________/
            | |                                    _________________| |_________________
            | |                                   |  _________________________________  |
4. Alignment (single-end)           5. Alignment 1 (paired-end)         6. Alignment 2 (paired-end)
            \ \                                    \ \_______________________________/ /
             \ \                                    \_______________   _______________/
              \ \                                                   | |                                    
               \ \                                            7. Merge (paired-end)
                \ \________________________________________________/ /
                 \______________________   _________________________/
                                        | |
                                8. Mark duplicates
            ____________________________| |
           |  __________________________  |
           | |                          | |
 Pre-filtering statistics          14. Filtering
 __________| |____________              | |
|                         |  Post-filtering statistics
|9. Extract unmapped reads|  ___________| |___________
|     10. idxstats        | |                         |
|     11. flagstat        | |      15. coverage       |
|     12. coverage        | |      16. idxstat        |
|      13. stats          | |      17. coverage       |
|_________________________| |   18. change in reads   |
                            |        19. stats        |
                            |      20. qualimap       |
                            |     21. read-depth      |
                            |___________   ___________|
                                        | |
                                    22. Clean-up

EOF
}


# ------------------ Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 's:d:a:r:ufj:t:hv' OPTION; do
    case "$OPTION" in
        s)
            if [ -d "$OPTARG" ]; then
                SD="$(readlink -f "$OPTARG")"
            else
                echo -e "\n$OPTARG does not seem to be a directory\n\nIf unsure of how to proceed run: $(basename "$script") -h\n"
                exit 1
            fi
            ;;
        d)
            if [ -d "$OPTARG" ]; then
                WD="$(readlink -f "$OPTARG")"
            else
                echo -e "\n$OPTARG does not seem to be a directory\n\nIf unsure of how to proceed run: $(basename "$script") -h\n"
                exit 1
            fi
            ;;
        a)
            if [ "$OPTARG" == "mem" ] || [ "$OPTARG" == "aln" ]; then
                algo="$OPTARG"
            else
                echo -e "\n-a must be either 'mem' [default] or 'aln'\n\nIf unsure of how to proceed run: $(basename "$script") -h\n"
				exit 1
	    	fi
	    	;;
        r)
            if [[ "$OPTARG" == *.@(fna|fa|fasta) ]]; then
                RG=$(readlink -f "$OPTARG")
            else
                echo -e "\nSupplied reference genome must be in one of the following formats .fna, .fa or .fasta\n"
                exit 1
            fi
            ;;
        u)
            single_sample="Y"
            ;;
        f)
            force_overwrite="Y"
            ;;
        j)
            if [ "$OPTARG" == "list" ]; then
                jobgraph
                exit 0
            else
                for (( i=1; i<="$num_of_jobs"; i++ )); do
                    joblist["$i"]=false
                done
                if [ ${#OPTARG} -gt 2 ]; then
                    min=${OPTARG%-*}
                    max=${OPTARG#*-}
                    for (( i="$min"; i<="$max"; i++ )); do
                        joblist["$i"]=true
                    done
                else
                    joblist["$OPTARG"]=true
                fi
            fi
            ;;
        t)
            timefactor="$OPTARG"
            ;;
        h)
            usage
            exit 0
            ;;
        v)
            version
            exit 0
            ;;
        ?)
            usage
            exit 2
            ;;
    esac
done

# Removes "/" from the end of directory path
if [ "${SD: -1}" == "/" ]; then
    length=${#SD}
    SD=${SD::length - 1}
fi

# Removes "/" from the end of directory path
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

# Function to flexibly control the dependencies of each job
dependency(){
    depend_on=
    for id in $1; do
         if [ -n "${jobidlist[$id]}" ]; then
            if [ -z "$depend_on" ]; then
                depend_on="${jobidlist[$id]}"
            else
                depend_on="$depend_on:${jobidlist[$id]}"
            fi
        fi
    done
    [ -z "$depend_on" ] && depend_on="-1"
    echo "$depend_on"
}

# Function used to adjust time for jobs
timer(){
    adjbase=$(awk -v adjustment="$1" -v base="$2" -v timefactor="$timefactor" 'BEGIN { print int( ( base * adjustment ) * 1.50 * timefactor ) }')
    if [ $((adjbase)) -lt "$3" ]; then 
        awk -v adjbase="$adjbase" -v static="$3" 'BEGIN { print int( adjbase + static ) }'
    elif [ $((adjbase)) -ge 10080 ]; then
        echo "10080"
    else
        echo "$adjbase"
    fi
}

# 01 AdapterRemoval jobs
job1(){
    echo "$(basename "$sample") is sent to the queue as a single-end$age sample" >> "$logfile"
    jobidlist[1]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 720 120)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "")" \
                --output="$out"/01_single_adapterremoval-%j.out \
                "$script_path"/modules/02_data_prep_01_single_adapterremoval.sh "$WD" "$sample" "$temp")
    echo -e "\t'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[1]}" >> "$logfile"
}

job23(){
    echo "$(basename "$sample") is sent to the queue as a single-end$age sample" >> "$logfile"
    jobidlist[23]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 720 120)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "")" \
                --output="$out"/01_singleseveral_adapterremoval-%j.out \
                "$script_path"/modules/02_data_prep_01_singleseveral_adapterremoval.sh "$WD" "$sample" "$temp")
    echo -e "\t'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[23]}" >> "$logfile"
}

job2(){
    echo "$(basename "$sample") is sent to the queue as a pair-end$age sample" >> "$logfile"
    jobidlist[2]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 720 120)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "")" \
                --output="$out"/01_paired_adapterremoval-%j.out \
                "$script_path"/modules/02_data_prep_01_paired_adapterremoval.sh "$WD" "$sample" "$temp")
    echo -e "\t'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[2]}" >> "$logfile"
}

job3(){
    echo "$(basename "$sample") is sent to the queue as a pair-end$age sample" >> "$logfile"
    jobidlist[3]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 820 180)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "")" \
                --output="$out"/01_pairedseveral_adapterremoval-%j.out \
                "$script_path"/modules/02_data_prep_01_pairedseveral_adapterremoval.sh "$WD" "$sample" "$temp")
    echo -e "\t'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[3]}" >> "$logfile"
}

# 02 Alignemt Jobs
job4(){
    # Base time increased from 1800 to 3400 to accomodate increase time need for historic single-end sample
    jobidlist[4]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 3400 120)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "1 23")" \
                --output="$out"/02_single_alignment-%j.out \
                "$script_path"/modules/02_data_prep_02_single_alignment.sh "$RG" "$WD" "$sample" "$temp" "$algo" "$age")
    echo -e "\t'Alignment' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[4]}" >> "$logfile"
}

job5(){
    jobidlist[5]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 2000 180)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "2 3")" \
                --output="$out"/02_01_paired_alignment-%j.out \
                "$script_path"/modules/02_data_prep_02_01_paired_alignment.sh "$RG" "$WD" "$sample" "$temp" "$algo" "$age")
    echo -e "\t'Alignment of paired ends' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[5]}" >> "$logfile"
}

job6(){
    jobidlist[6]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 2000 180)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "2 3")" \
                --output="$out"/02_02_paired_alignment-%j.out \
                "$script_path"/modules/02_data_prep_02_02_paired_alignment.sh "$RG" "$WD" "$sample" "$temp" "$algo" "$age")
    echo -e "\t'Alignment of collapsed single end' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[6]}" >> "$logfile"
}

# 03 Merge paired alignment
job7(){
    jobidlist[7]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 1200 180)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "5 6")" \
                --output="$out"/03_paired_merge-%j.out \
                "$script_path"/modules/02_data_prep_03_paired_merge.sh "$WD" "$sample")
    echo -e "\t'Merging' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[7]}" >> "$logfile"
}

# 04 Mark duplicates
job8(){
    jobidlist[8]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 1800 120)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "4 7")" \
                --output="$out"/04_markduplicates-%j.out \
                "$script_path"/modules/02_data_prep_04_markduplicates.sh "$WD" "$sample" "$temp")
    echo -e "\t'Marking duplicates' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[8]}" >> "$logfile"
}

# 05 Extract unmapped reads
job9(){
    jobidlist[9]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 120 120)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "8")" \
                --output="$out"/05_unmapped-%j.out \
                "$script_path"/modules/02_data_prep_05_unmapped.sh "$SD" "$WD" "$sample")
    echo -e "\t'Extract unmapped' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[9]}" >> "$logfile"
}

# 06 Pre-filtering statistics
job10(){
    jobidlist[10]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 100 100)" \
                --mem=48G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "8")" \
                --output="$out"/06_01_prestats-%j.out \
                "$script_path"/modules/02_data_prep_06_01_prestats.sh "$WD" "$sample")
    echo -e "\t'Pre-filtering statistics (idxstats)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[10]}" >> "$logfile"
}

job11(){
    jobidlist[11]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 100 100)" \
                --mem=48G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "8")" \
                --output="$out"/06_02_prestats-%j.out \
                "$script_path"/modules/02_data_prep_06_02_prestats.sh "$WD" "$sample")
    echo -e "\t'Pre-filtering statistics (flagstat)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[11]}" >> "$logfile"
}

job12(){
    jobidlist[12]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 420 180)" \
                --mem=20G \
                --cpus-per-task=2 \
                --dependency=afterany:"$(dependency "8")" \
                --output="$out"/06_03_prestats-%j.out \
                "$script_path"/modules/02_data_prep_06_03_prestats.sh "$WD" "$sample")
    echo -e "\t'Pre-filtering statistics (coverage)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[12]}" >> "$logfile"
}

job13(){
    jobidlist[13]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 180 180)" \
                --mem=48G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "8")" \
                --output="$out"/06_04_prestats-%j.out \
                "$script_path"/modules/02_data_prep_06_04_prestats.sh "$WD" "$sample")
    echo -e "\t'Pre-filtering statistics (stats)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[13]}" >> "$logfile"
}

# 07 Filtering
job14(){
    jobidlist[14]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 1560 220)" \
                --mem=64G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "8")" \
                --output="$out"/07_filtering-%j.out \
                "$script_path"/modules/02_data_prep_07_filtering.sh "$SD" "$WD" "$sample")
    echo -e "\t'Filtering' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[14]}" >> "$logfile"
}

# 08 Post-filtering statistics
job15(){
    jobidlist[15]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 120 120)" \
                --mem=48G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "14")" \
                --output="$out"/08_01_poststats-%j.out \
                "$script_path"/modules/02_data_prep_08_01_poststats.sh "$SD" "$WD" "$sample")
    echo -e "\t'Post-filtering statistics (flagstat)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[15]}" >> "$logfile"
}

job16(){
    jobidlist[16]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 120 120)" \
                --mem=48G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "14")" \
                --output="$out"/08_02_poststats-%j.out \
                "$script_path"/modules/02_data_prep_08_02_poststats.sh "$SD" "$WD" "$sample")
    echo -e "\t'Post-filtering statistics (idxstats)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[16]}" >> "$logfile"
}

job17(){
    jobidlist[17]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 420 180)" \
                --mem=20G \
                --cpus-per-task=2 \
                --dependency=afterany:"$(dependency "14")" \
                --output="$out"/08_03_poststats-%j.out \
                "$script_path"/modules/02_data_prep_08_03_poststats.sh "$SD" "$WD" "$sample")
    echo -e "\t'Post-filtering statistics (coverage)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[17]}" >> "$logfile"
}

job18(){
    jobidlist[18]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 420 180)" \
                --mem=48G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "14")" \
                --output="$out"/08_04_poststats-%j.out \
                "$script_path"/modules/02_data_prep_08_04_poststats.sh "$SD" "$WD" "$sample")
    echo -e "\t'Post-filtering statistics (readchange)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[18]}" >> "$logfile"
}

job19(){
    jobidlist[19]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 180 180)" \
                --mem=48G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "14")" \
                --output="$out"/08_05_poststats-%j.out \
                "$script_path"/modules/02_data_prep_08_05_poststats.sh "$SD" "$WD" "$sample")
    echo -e "\t'Post-filtering statistics (stats)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[19]}" >> "$logfile"
}

job20(){
    jobidlist[20]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 300 180)" \
                --mem=160G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "14")" \
                --output="$out"/08_06_poststats-%j.out \
                "$script_path"/modules/02_data_prep_08_06_qualimap.sh "$RD" "$SD" "$sample")
    echo -e "\t'Post-filtering statistics (qualimap)' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[20]}" >> "$logfile"
}

job21(){
    jobidlist[21]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time="$(timer "$adjustment" 60 60)" \
                --mem=200G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "14")" \
                --output="$out"/08_07_poststats-%j.out \
                "$script_path"/modules/02_data_prep_08_07_depth.sh "$SD" "$WD" "$sample" "$script_path")
    echo -e "\t'Post-filtering read depth' job has been submitted for $(basename "$sample") -- Job ID: ${jobidlist[21]}\n" >> "$logfile"
}

# 09 Clean up
job22(){
    jobidlist[22]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time=30 \
                --mem=2G \
                --cpus-per-task=1 \
                --dependency=afterany:"$(dependency "15 16 17 18 19 20 21")" \
                --output=/dev/null \
                --error=/dev/null \
                "$script_path"/modules/02_data_prep_09_cleanup.sh "$sample" "$out" "$logfile" "$starttime")
}

# Function for sample processing
sample_processing(){
	# Checks if sample folder is empty
	if [ "$(ls -A "$sample")" ]; then	
		# Checks whether a .bam file already exists within sample folder, indicating samples have already been processed        
		for bam in "$sample"/*filtered.bam; do
			if [ ! -e "$bam" ] || [ "$force_overwrite" == "Y" ]; then
				# Checks sample file size
				filesize=
                counter=0
				for fna in "$sample"/*.fq.gz; do
                    ((counter++))
					if [ ! "$filesize" ]; then
                        filesize=$(wc -c < "$fna")
                    fi
				done
                if [ $((counter)) -gt 2 ]; then
                    filesize=$(( filesize * (counter / 2) ))
                fi
				# Change requested time on nodes depending on filesize comparative to R1 Ocin_NYS-F
				adjustment=$(awk -v filesize="$filesize" 'BEGIN { print ( filesize / 93635424798 ) }')
				# Creates sample directory in species directory if none exist
				sampledir="$WD"/"$(basename "$sample")"
				[[ -d "$sampledir" ]] || mkdir -m 775 "$sampledir"
				# Creates pre- and post-filtering directory in sample directory if none exist
				[[ -d "$sampledir"/pre_filter_stats ]] || mkdir -m 775 "$sampledir"/pre_filter_stats
				[[ -d "$sampledir"/post_filter_stats ]] || mkdir -m 775 "$sampledir"/post_filter_stats
				# Create folder for STDOUT files generated sbatch
				out="$sampledir"/out
				[[ -d "$out" ]] || mkdir -m 775 "$out"
				# Clean out old out files
				[ ! "$(ls -A "$out")" ] || rm "$out"/*
				# Checks if currently working with museomics samples
				if [[ "$SD" == *"museomics"* ]]; then
					# Checks if sample directory is pre-2000 (historic) or post-2000 (modern) and sets the algrorithm parameter correspondingly
					if [ $((${sample: -4})) -gt 2000 ]; then
						# Sample is modern
						age=" modern"
						algo="mem"
						queue
					else
						# Sample is historic
						age=" historic"
						algo="aln"
						queue
					fi
				else
					# Sample is custom
					age=""
					queue
				fi
			else
				echo -e "$(basename "$sample") already contains a .bam file, $(basename "$bam"), and is skipped\n" >> "$logfile"
			fi
		done
	else
		echo -e "$(basename "$sample") is an empty directory and is skipped\n" >> "$logfile"
	fi
}

# Function used to queue jobs
queue(){
	for _2 in "$sample"/*_2.fq.gz; do
        if [ ! -e "$_2" ]; then
            count=$(find -L "$sample"/ -maxdepth 1 -type f -name '*_1.fq.gz' | wc -l)
            # Single-end
            if [ $((count)) -eq 1 ]; then
                # AdapterRemoval
                # Executes job if indicated in joblist
                [ "${joblist[1]}" = true ] && job1
            # Several single end fiels
            elif [ $((count)) -gt 1 ]; then
                # AdapterRemoval
                # Executes job if indicated in joblist
                [ "${joblist[23]}" = true ] && job23
            fi
            # Aligning to reference
            # Executes job if indicated in joblist
            [ "${joblist[4]}" = true ] && job4

            # Marking duplicates
            # Executes job if indicated in joblist
            [ "${joblist[8]}" = true ] && job8
        else
            count=$(find -L "$sample"/ -maxdepth 1 -type f -name '*_2.fq.gz' | wc -l)
            # Paired-end
            if [ $((count)) -eq 1 ]; then
                # AdapterRemoval
                # Executes job if indicated in joblist
                [ "${joblist[2]}" = true ] && job2
            # Several paired end fiels
            elif [ $((count)) -gt 1 ]; then
                #AdapterRemoval
                # Executes job if indicated in joblist
                [ "${joblist[3]}" = true ] && job3
            fi
            # Aligning to reference
            # Executes job if indicated in joblist
            [ "${joblist[5]}" = true ] && job5
            # Executes job if indicated in joblist
            [ "${joblist[6]}" = true ] && job6

            # Merging of alignment files
            # Executes job if indicated in joblist
            [ "${joblist[7]}" = true ] && job7

            # Marking duplicates
            # Executes job if indicated in joblist
            [ "${joblist[8]}" = true ] && job8
        fi
        break
    done
	# Extract unmapped reads
	# Executes job if indicated in joblist
    [ "${joblist[9]}" = true ] && job9

	# Statistics pre-filtering
	# Executes job if indicated in joblist
    [ "${joblist[10]}" = true ] && job10
    # Executes job if indicated in joblist
    [ "${joblist[11]}" = true ] && job11
    # Executes job if indicated in joblist
    [ "${joblist[12]}" = true ] && job12
    # Executes job if indicated in joblist
    [ "${joblist[13]}" = true ] && job13

	# Removal of duplicates, unmapped reads and low quality mappings
	# Executes job if indicated in joblist
    [ "${joblist[14]}" = true ] && job14

	# Statistics post-filtering
	# Executes job if indicated in joblist
    [ "${joblist[15]}" = true ] && job15
    # Executes job if indicated in joblist
    [ "${joblist[16]}" = true ] && job16
    # Executes job if indicated in joblist
    [ "${joblist[17]}" = true ] && job17
    # Executes job if indicated in joblist
    [ "${joblist[18]}" = true ] && job18
    # Executes job if indicated in joblist
    [ "${joblist[19]}" = true ] && job19
    # Executes job if indicated in joblist
    [ "${joblist[20]}" = true ] && job20
    # Executes job if indicated in joblist
    [ "${joblist[21]}" = true ] && job21

	# Clean up of empty stdout files
	# Executes job if indicated in joblist
    [ "${joblist[22]}" = true ] && job22
}

# ----------------- Main -------------------------------------------------

# Finds references genome for designated species or exits if it cannot be found
if [ ! "$RG" ]; then
    for ref in /faststorage/project/EcoGenetics/BACKUP/reference_genomes/"$(basename "$SD")"/*.@(fna|fa|fasta); do
        if [ ! -e "$ref" ]; then
            echo -e "\nCannot locate reference genome $ref, format: .fna, .fa, .fasta\n"
            exit 1
        else
            RG=$ref
            break
        fi
    done
fi

# Exit with error message if reference genome is not indexed
for index in "$(dirname "$RG")"/*.ann; do
    if [ ! -e "$index" ]; then
		echo -e "\nDesignated reference genome does not seem to be indexed, $RG\n"
		exit 1
    fi
done

# Creates working directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"
if [ ! -w "$WD" ]; then
    echo -e "\nYou do not have write permission in the selected working directory: $WD\n"
    exit 3
fi

# Creates temp directory in working directory if none exist
temp="$WD"/temp
[ -d "$temp" ] || mkdir -m 775 "$temp"

# Creates data directory in working directory if none exist
WD="$WD"/"$datadir"
[ -d "$WD" ] || mkdir -m 775 "$WD"

# Creates species abbreviation from species name
species_name="$(basename "$SD")"
genus=${species_name%_*}; genus=${genus::3}; genus=${genus^}
species=${species_name#*_}; species=${species::3}; species=${species^}
speciesabbr="$genus""$species"

# Creates species directory in data_preparation directory if none exist
WD="$WD"/"$species_name"
[ -d "$WD" ] || mkdir -m 775 "$WD"

# Removes old log files and creates new log file
if [[ -e "$WD"/log_"$speciesabbr"_data_preparation.txt ]]; then
    rm "$WD"/log_"$speciesabbr"_data_preparation.txt
fi
logfile="$WD"/log_"$speciesabbr"_data_preparation.txt
touch "$logfile"
starttime="$(date +"%Y-%m-%d")"
echo -e "----- $starttime -----\n" >> "$logfile"

# Makes sure all associated modules are executable
for script in "$script_path"/modules/*; do
	[ ! -x "$script" ] && chmod -R u+rwx "$script_path"/modules && break
done

# Checks whether or not to run on single sample
if [ "$single_sample" == "N" ]; then
	# Multiple samples
	# Loops through all sample folders within species specific sample directory
	for sample in "$SD"/*; do
        if [ ! -d "$sample" ]; then
            # Ignores anything in is species directory which is not itself a directory
            continue
        else
            # Makes sure list containing job IDs is empty
            jobidlist=()
            sample_processing
        fi
	done
elif [ "$single_sample" == "Y" ]; then
	# Single sample
    # Makes sure list containing job IDs is empty
    jobidlist=()
	sample_processing
fi

echo -e "\nLog file can be found here: $logfile\n"

exit 0