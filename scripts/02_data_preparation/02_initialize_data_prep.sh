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

# ----------------- Description ------------------------------------------

# Script for initializing data preparation procedure for sequence data
# Tested and working using:
# EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary

fi

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 02_initialize_data_prep.sh [PARAMETERS] [OPTIONS]

This script is used for initializing the standardized data preparation procedure
for sequence data at the Center for EcoGenetics.
Can be run on the frontend, as its resource-demands are fairly, or run in 
conjunction with 'sbatch'or 'srun'.

If the specified species sample directory is within the 'museomics' directory,
the script will automatically try to detect whether the sample directory is
contemporary or historical and automatically apply the corresponding algorithm.
This means that if your samples are in the museomics directory you do not need
to specify an algorithm.

PARAMETERS (must be assigned):
    -s  DIRECTORY       Species specific sample directory

OPTIONS:
    -d  DIRECTORY       Working directory. If not assigned, will use current 
                        working directory [default]
    -a  ALGORITHM       Choice of algorithm to be used during alignment.
                        'mem' (>70MB, contemporary samples)[default] or 
                        'aln' (<70MB, historic samples)
    -m  INTEGER         Amount of memory to be used by each CPU. 8 [default]
    -c  INTEGER         Number of CPUs to be used. 8 [default]
    -u                  Run only on a single sample. -s then needs to lead to
                        the specific sample
    -f                  Force run even if target sample directory contains .bam
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

# Working directory
WD="$(readlink -f "$PWD")"

# Algorithm to use during alignment: mem (>70MB, contemporary samples) [default] or aln (<70MB, historic samples)
algo="mem"

# Define memory per cpu in G (must be integer)
memory="8"

# Define number of cpus to be used (must be integer)
cpus="8"

# Option to only run a single sample
single_sample="N"

# Option to force run even if target sample directory contains .bam
force_overwrite="N"

# Name of created data directory under the working directory
dataprep="02_data_preparation"

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

# ----------------- Script Flag Processing -------------------------------

while getopts 's:d:a:m:c:ufh' OPTION; do
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
	    	if [[ $OPTARG == "mem" ]] || [[ $OPTARG == "aln" ]]; then
				algo="$OPTARG"
	    	else
				echo -e "\nERROR: -a must be either 'mem' [default] or 'aln'\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
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
		
		# Checks whether a .bam file already exists within sample folder, indicating samples have already been processed        
		for bam in "$sample"/*.bam; do
			if [ ! -e "$bam" ] || [ "$force_overwrite" == "Y" ]; then
			
				# Checks sample file size
				filesize=
				for fna in "$sample"/*.fq.gz; do
					[ ! "$filesize" ] && filesize=$(wc -c < "$fna") && break
				done

				# Change requested time on nodes depending on filesize comparative to R1 Ocin_NYS-F
				adjustment=$(awk -v filesize="$filesize" 'BEGIN { print ( filesize / 93635424798 + 0.1) }')

				# Creates sample directory in species directory if none exist
				sampledir="$speciesdir"/"$(basename "$sample")"
				[[ -d "$sampledir" ]] || mkdir -m 775 "$sampledir"

				# Creates pre- and post-filtering directory in sample directory if none exist
				[[ -d "$sampledir"/pre_filter_stats ]] || mkdir -m 775 "$sampledir"/pre_filter_stats
				[[ -d "$sampledir"/post_filter_stats ]] || mkdir -m 775 "$sampledir"/post_filter_stats

				# Create folder for STDOUT files generated sbatch
				stdoutput="$sampledir"/stdout_sbatch
				[[ -d "$stdoutput" ]] || mkdir -m 775 "$stdoutput"

				# Checks if currently working with museomics samples
				if [[ "$SD" == *"museomics"* ]]; then

					# Checks if sample directory is pre-2000 (historic) or post-2000 (modern) and sets the algrorithm parameter correspondingly
					if [ "$((${sample: -4}))" -gt 2000 ]; then

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

# Function used to adjust time for jobs
# timer()
# {
# 	awk -v adjustment="$1" -v base="$2" -v static="$3" 'BEGIN { print int( base * adjustment + static) }'
# }
timer()
{
    adjbase=$(awk -v adjustment="$1" -v base="$2" 'BEGIN { print int( base * adjustment ) }')
	
    if [ "$adjbase" -lt "$3" ]; then 
        
        awk -v adjbase="$adjbase" -v static="$3" 'BEGIN { print int( adjbase + static ) }'
    
    else

        echo "$adjbase"
    
    fi
}

# Function used to queue jobs
queue()
{
	# Check the number of .fq.gz files in sample directory (assumed to be indicative of whether sample is single- or paired-end)
	count=$(find "$sample"/ -maxdepth 1 -type f -name '*.fq.gz' | wc -l)
	if [ "$count" == 1 ]; then
			    
		# Single-end

		echo "$(basename "$sample") is sent to the queue as a single-end$age sample" >> "$logfile"

		# AdapterRemoval
		jid1=$(sbatch \
			--parsable \
			--time="$(timer "$adjustment" 720 120)" \
			--mem-per-cpu="$memory"G \
			--cpus-per-task="$cpus" \
			--output="$stdoutput"/"$(basename "$sample")"_01-%j.out \
			"$script_path"/modules/02_01_single_adapterremoval.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
			    
		echo -e "\t'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: $jid1" >> "$logfile"

		# Aligning to reference
		jid2=$(sbatch \
			--parsable \
			--time="$(timer "$adjustment" 1800 120)" \
			--mem-per-cpu="$memory"G \
			--cpus-per-task="$cpus" \
			--output="$stdoutput"/"$(basename "$sample")"_02-%j.out \
			--dependency=afterany:"$jid1" \
			"$script_path"/modules/02_02_single_alignment.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")

		echo -e "\t'Alignment' job has been submitted for $(basename "$sample") -- Job ID: $jid2" >> "$logfile"

		# Marking duplicates
		jid4=$(sbatch \
			--parsable \
			--time="$(timer "$adjustment" 1800 120)" \
			--mem-per-cpu="$memory"G \
			--cpus-per-task="$cpus" \
			--output="$stdoutput"/"$(basename "$sample")"_04-%j.out \
			--dependency=afterany:"$jid2" \
			"$script_path"/modules/02_04_markduplicates.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")

		echo -e "\t'Marking duplicates' job has been submitted for $(basename "$sample") -- Job ID: $jid4" >> "$logfile"
			
	else

		# Paired-end

		echo "$(basename "$sample") is sent to the queue as a pair-end$age sample" >> "$logfile"

		if [ "$count" == 2 ];then	    

			# AdapterRemoval
			jid1=$(sbatch \
				--parsable \
				--time="$(timer "$adjustment" 720 120)" \
				--mem-per-cpu="$memory"G \
				--cpus-per-task="$cpus" \
				--output="$stdoutput"/"$(basename "$sample")"_01-%j.out \
				"$script_path"/modules/02_01_paired_adapterremoval.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
					
			echo -e "\t'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: $jid1" >> "$logfile"
		
		else

			#AdapterRemoval
			jid1=$(sbatch \
				--parsable \
				--time="$(timer "$adjustment" 720 150)" \
				--mem-per-cpu="$memory"G \
				--cpus-per-task="$cpus" \
				--output="$stdoutput"/"$(basename "$sample")"_01-%j.out \
				"$script_path"/modules/02_01_pairedseveral_adapterremoval.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
					
			echo -e "\t'AdapterRemoval' job has been submitted for $(basename "$sample") -- Job ID: $jid1" >> "$logfile"

		fi

		# Aligning to reference
		jid2_1=$(sbatch \
			--parsable \
			--time="$(timer "$adjustment" 1800 120)" \
			--mem-per-cpu="$memory"G \
			--cpus-per-task="$cpus" \
			--output="$stdoutput"/"$(basename "$sample")"_02_01-%j.out \
			--dependency=afterany:"$jid1" \
			"$script_path"/modules/02_02_paired_01_alignment.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
		jid2_2=$(sbatch \
			--parsable \
			--time="$(timer "$adjustment" 1800 120)" \
			--mem-per-cpu="$memory"G \
			--cpus-per-task="$cpus" \
			--output="$stdoutput"/"$(basename "$sample")"_02_02-%j.out \
			--dependency=afterany:"$jid1" \
			"$script_path"/modules/02_02_paired_02_alignment.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")

		echo -e "\t'Alignment of paired ends' job has been submitted for $(basename "$sample") -- Job ID: $jid2_1" >> "$logfile"
		echo -e "\t'Alignment of collapsed single end' job has been submitted for $(basename "$sample") -- Job ID: $jid2_2" >> "$logfile"

		# Merging of alignment files
		jid3=$(sbatch \
			--parsable \
			--time="$(timer "$adjustment" 1200 120)" \
			--mem-per-cpu="$memory"G \
			--cpus-per-task="$cpus" \
			--output="$stdoutput"/"$(basename "$sample")"_03-%j.out \
			--dependency=afterany:"$jid2_1":"$jid2_2" \
			"$script_path"/modules/02_03_paired_merge.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")

		echo -e "\t'Merging' job has been submitted for $(basename "$sample") -- Job ID: $jid3" >> "$logfile"

		# Marking duplicates
		jid4=$(sbatch \
			--parsable \
			--time="$(timer "$adjustment" 1800 120)" \
			--mem-per-cpu="$memory"G \
			--cpus-per-task="$cpus" \
			--output="$stdoutput"/"$(basename "$sample")"_04-%j.out \
			--dependency=afterany:"$jid3" \
			"$script_path"/modules/02_04_markduplicates.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
			    
		echo -e "\t'Marking duplicates' job has been submitted for $(basename "$sample") -- Job ID: $jid4" >> "$logfile"

	fi
			
	# Statistics pre-filtering
	jid5_1=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 60 60)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_05_01-%j.out \
		--dependency=afterany:"$jid4" \
		"$script_path"/modules/02_05_01_prestats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
	jid5_2=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 60 60)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_05_02-%j.out \
		--dependency=afterany:"$jid4" \
		"$script_path"/modules/02_05_02_prestats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
	jid5_3=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 360 120)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_05_03-%j.out \
		--dependency=afterany:"$jid4" \
		"$script_path"/modules/02_05_03_prestats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
	jid5_4=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 60 60)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_05_04-%j.out \
		--dependency=afterany:"$jid4" \
		"$script_path"/modules/02_05_04_prestats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
			
	echo -e "\t'Pre-filtering statistics (idxstats)' job has been submitted for $(basename "$sample") -- Job ID: $jid5_1" >> "$logfile"
	echo -e "\t'Pre-filtering statistics (flagstat)' job has been submitted for $(basename "$sample") -- Job ID: $jid5_2" >> "$logfile"
	echo -e "\t'Pre-filtering statistics (coverage)' job has been submitted for $(basename "$sample") -- Job ID: $jid5_3" >> "$logfile"
	echo -e "\t'Pre-filtering statistics (stats)' job has been submitted for $(basename "$sample") -- Job ID: $jid5_4" >> "$logfile"


	# Removal of duplicates, unmapped reads and low quality mappings
	jid6=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 1440 120)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_06-%j.out \
		--dependency=afterany:"$jid5_1":"$jid5_2":"$jid5_3":"$jid5_4" \
		"$script_path"/modules/02_06_filtering.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")

	echo -e "\t'Filtering' job has been submitted for $(basename "$sample") -- Job ID: $jid6" >> "$logfile"

	# Statistics post-filtering
	jid7_1=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 60 60)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_07_01-%j.out \
		--dependency=afterany:"$jid6" \
		"$script_path"/modules/02_07_01_poststats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
	jid7_2=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 60 60)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_07_02-%j.out \
		--dependency=afterany:"$jid6" \
		"$script_path"/modules/02_07_02_poststats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
	jid7_3=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 360 120)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_07_03-%j.out \
		--dependency=afterany:"$jid6" \
		"$script_path"/modules/02_07_03_poststats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
	jid7_4=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 360 120)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_07_04-%j.out \
		--dependency=afterany:"$jid6" \
		"$script_path"/modules/02_07_04_poststats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
	jid7_5=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 60 60)" \
		--mem-per-cpu="$memory"G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_07_05-%j.out \
		--dependency=afterany:"$jid6" \
		"$script_path"/modules/02_07_05_poststats.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
	jid7_6=$(sbatch \
		--parsable \
		--time="$(timer "$adjustment" 240 120)" \
		--mem-per-cpu=20G \
		--cpus-per-task="$cpus" \
		--output="$stdoutput"/"$(basename "$sample")"_07_06-%j.out \
		--dependency=afterany:"$jid6" \
		"$script_path"/modules/02_07_06_qualimap.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")

	echo -e "\t'Post-filtering statistics (flagstat)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_1" >> "$logfile"
	echo -e "\t'Post-filtering statistics (idxstats)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_2" >> "$logfile"
	echo -e "\t'Post-filtering statistics (coverage)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_3" >> "$logfile"
	echo -e "\t'Post-filtering statistics (readchange)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_4" >> "$logfile"
	echo -e "\t'Post-filtering statistics (stats)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_5" >> "$logfile"
	echo -e "\t'Post-filtering statistics (qualimap)' job has been submitted for $(basename "$sample") -- Job ID: $jid7_6\n" >> "$logfile"

	# Clean up of empty stdout files
	jid8=$(sbatch \
		--output=/dev/null \
		--error=/dev/null \
		--dependency=afterany:"$jid7_1":"$jid7_2":"$jid7_3":"$jid7_4":"$jid7_5":"$jid7_6" \
		"$script_path"/modules/02_08_cleanup.sh "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo")
}

# ----------------- Script Queue -----------------------------------------

# Holds path to script directory
script_path=$(dirname "$script_path")

# Finds references genome for designated species or exits if it cannot be found
shopt -s nullglob extglob
for ref in /faststorage/project/EcoGenetics/BACKUP/reference_genomes/"$(basename "$SD")"/*.@(fna|fa|fasta); do
    if [ ! -e "$ref" ]; then
		echo -e "\nCannot locate reference genome $ref, format: .fna, .fa, .fasta\n"
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
