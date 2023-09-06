#!/bin/bash

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

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 03_initialize.sh [PARAMETERS] [OPTIONS]

This script is used for initializing the calculation of Site Frequency Spectra
from sequence data at the Center for EcoGenetics.
Can be run on the frontend, as its resource-demands are fairly, or run in 
conjunction with 'sbatch'or 'srun'.

PARAMETERS (must be assigned):
    -s  FILE            Sample .mpileup file
    -a  INTEGER         Number of alleles. For diploid organisms 
                        = 2*[number of individuals]

OPTIONS:
    -d  DIRECTORY       Working directory. If not assigned, will use current 
                        working directory [default]
    -f                  Force run even if target file already has a matching 
                        .spectrum file
    -t                  Keep created temporary files. Mostly for testing
                        purposes
    -h                  Show this message


EOF
}

# ----------------- Configuration ----------------------------------------

# Sample .mpileup file, abosolute path
S=

# Number of alleles
A=

# Working directory
WD="$(readlink -f "$PWD")"

# Option to force run even if target file already has a matching .spectrum file
force_overwrite="N"

# Option to keep all created temp files
keep_temp="N"

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

# ----------------- Script Flag Processing -------------------------------

while getopts 's:a:d:fth' OPTION; do
    case "$OPTION" in
		s)
	    	if [ -f "$OPTARG" ]; then
				S="$(readlink -f "$OPTARG")"
	    	else
				echo -e "\nERROR: $OPTARG does not seem to be a file\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
				exit 1
	    	fi
	    	;;
        a)
	    	if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
				A="$OPTARG"
	    	else
				echo -e "\nERROR: $OPTARG is not an integer\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
				exit 1
	    	fi
	    	;;
		d)
	    	WD="$(readlink -f "$OPTARG")"
	    	;;
		f)
	    	force_overwrite="Y"
	    	;;
		t)
	    	keep_temp="Y"
	    	;;
		h)
	    	usage
	    	exit 0
	    	;;
		?)
	    	usage
	    	exit 1
	    	;;
    esac
done

# Removes "/" from the end of path for working directory
if [ "${WD: -1}" == "/" ]; then
    length=${#WD}
    WD=${WD::length - 1}
fi

# Tests that sample .mpileup file has been supplied
if [[ -z $S ]]; then
    usage
    exit 1
fi

# Tests that allele number has been supplied
if [[ -z $A ]]; then
    usage
    exit 1
fi

# ----------------- Functions --------------------------------------------

# Function used to adjust time for jobs
timer()
{
    adjbase=$(awk -v adjustment="$1" -v base="$2" 'BEGIN { print int( base * adjustment ) }')
	
    if [ "$adjbase" -lt "$3" ]; then 
        awk -v adjbase="$adjbase" -v static="$3" 'BEGIN { print int( adjbase + static ) }'
    else
        echo "$adjbase"
    fi
}

filechecker()
{
    runlist=()

    # Checks for .mpileup parts
    for mpileup_part in "$temp"/"$samplename"_*.mpileup; do
        if [ -e "$mpileup_part" ]; then
            runlist+=("1")
            break
        else
            runlist+=("0")
            break
        fi
    done

    # Checks for filtered .mpileup parts
    for filtered_mpileup_part in "$temp"/"$samplename"_filtered_*.mpileup; do
        if [ -e "$filtered_mpileup_part" ]; then
            runlist+=("1")
            break
        else
            runlist+=("0")
            break
        fi
    done

    # Checks for .sync parts
    for sync_part in "$temp"/"$samplename"_*.sync; do
        if [ -e "$sync_part" ]; then
            runlist+=("1")
            break
        else
            runlist+=("0")
            break
        fi
    done

    # Checks for filtered .sync parts
    for filtered_sync_part in "$temp"/"$samplename"_filtered_*.sync; do
        if [ -e "$filtered_sync_part" ]; then
            runlist+=("1")
            break
        else
            runlist+=("0")
            break
        fi
    done

    # Checks for .spectrum parts
    for spectrum_part in "$temp"/"$samplename"_*.spectrum; do
        if [ -e "$spectrum_part" ]; then
            runlist+=("1")
            break
        else
            runlist+=("0")
            break
        fi
    done
}

job1()
{
    if [ "${runlist[0]}" == "0" ]; then
        jid1=$(sbatch \
            --parsable \
            --array=1-"$lines"%300 \
            --time="$(timer "$adjustment" 360 60)" \
            --mem-per-cpu=6G \
            --cpus-per-task=6 \
            --output="$out"/split_pileup-%a-%j-%A.out \
            "$script_path"/03_split_pileup.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
    else
        jid1=
    fi
}

job2()
{
    if [ "${runlist[1]}" == "0" ]; then
        if [ -n "$jid1" ];then
            jid2=$(sbatch \
                --parsable \
                --array=1-"$lines"%200 \
                --time="$(timer "$adjustment" 180 60)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --dependency=aftercorr:"$jid1" \
                --output="$out"/filter_coverage-%a-%j-%A.out \
                "$script_path"/03_filter_coverage.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        else
            jid2=$(sbatch \
                --parsable \
                --array=1-"$lines"%200 \
                --time="$(timer "$adjustment" 180 60)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --output="$out"/filter_coverage-%a-%j-%A.out \
                "$script_path"/03_filter_coverage.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        fi
    else
        jid2=
    fi
}

job3()
{
    if [ "${runlist[2]}" == "0" ]; then
        if [ -n "$jid2" ];then
            jid3=$(sbatch \
                --parsable \
                --array=1-"$lines"%300 \
                --time="$(timer "$adjustment" 60 30)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --dependency=aftercorr:"$jid2" \
                --output="$out"/mpileup_to_sync-%a-%j-%A.out \
                "$script_path"/03_mpileup_to_sync.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        else
            jid3=$(sbatch \
                --parsable \
                --array=1-"$lines"%300 \
                --time="$(timer "$adjustment" 60 30)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --output="$out"/mpileup_to_sync-%a-%j-%A.out \
                "$script_path"/03_mpileup_to_sync.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        fi
    else
        jid3=
    fi
}

job4()
{
    if [ "${runlist[3]}" == "0" ]; then
        if [ -n "$jid3" ];then
            jid4=$(sbatch \
                --parsable \
                --array=1-"$lines"%300 \
                --time="$(timer "$adjustment" 180 120)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --dependency=afterany:"$jid3" \
                --output="$out"/readsync-%a-%j-%A.out \
                "$script_path"/03_readsync.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        else
            jid4=$(sbatch \
                --parsable \
                --array=1-"$lines"%300 \
                --time="$(timer "$adjustment" 180 120)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --output="$out"/readsync-%a-%j-%A.out \
                "$script_path"/03_readsync.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        fi
    else
        jid4=
    fi
}

job5()
{
    if [ "${runlist[4]}" == "0" ]; then
        if [ -n "$jid4" ];then
            jid5=$(sbatch \
                --parsable \
                --array=1-"$lines"%100 \
                --time="$(timer "$adjustment" 180 180)" \
                --mem-per-cpu=15G \
                --cpus-per-task=4 \
                --dependency=afterany:"$jid4" \
                --output="$out"/calculate_sfs-%a-%j-%A.out \
                "$script_path"/03_calculate_sfs.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        else
            jid5=$(sbatch \
                --parsable \
                --array=1-"$lines"%100 \
                --time="$(timer "$adjustment" 180 180)" \
                --mem-per-cpu=15G \
                --cpus-per-task=4 \
                --output="$out"/calculate_sfs-%a-%j-%A.out\
                "$script_path"/03_calculate_sfs.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        fi
    else
        jid5=
    fi
}

job6()
{
    if [ "${runlist[4]}" == "0" ]; then
        if [ -n "$jid5" ];then
            jid6=$(sbatch \
                --parsable \
                --time="$(timer "$adjustment" 180 60)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --dependency=afterany:"$jid5" \
                --output="$out"/concat_sync-%j.out \
                "$script_path"/03_concat_sync.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
            jid7=$(sbatch \
                --parsable \
                --time="$(timer "$adjustment" 180 60)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --dependency=afterany:"$jid5" \
                --output="$out"/concat_filtered_sync-%j.out \
                "$script_path"/03_concat_filtered_sync.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
            jid8=$(sbatch \
                --parsable \
                --time="$(timer "$adjustment" 180 60)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --dependency=afterany:"$jid5" \
                --output="$out"/concat_spectrum-%j.out \
                "$script_path"/03_concat_spectrum.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        else
            jid6=$(sbatch \
                --parsable \
                --time="$(timer "$adjustment" 180 60)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --output="$out"/concat_sync-%j.out \
                "$script_path"/03_concat_sync.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
            jid7=$(sbatch \
                --parsable \
                --time="$(timer "$adjustment" 180 60)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --output="$out"/concat_filtered_sync-%j.out \
                "$script_path"/03_concat_filtered_sync.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
            jid8=$(sbatch \
                --parsable \
                --time="$(timer "$adjustment" 180 60)" \
                --mem-per-cpu=10G \
                --cpus-per-task=2 \
                --output="$out"/concat_spectrum-%j.out \
                "$script_path"/03_concat_spectrum.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8")
        fi
    else
        jid6=
        jid7=
        jid8=
    fi
}

# ----------------- Script Queue -----------------------------------------

# Establishes name base based on .mpileup file
mpileup=$(basename "$S") 
namebase=${mpileup%.*}

# Checks if spectrum file already exists and 
if [ "$force_overwrite" == "N" ] && [ -e "${S%.*}".spectrum ]; then

    echo -e "\nERROR: '$namebase.spectrum' seems to already exist. Remove file or use overwrite option.\n\nIf unsure of how to proceed run: $(basename "$script_path") -h\n"
	exit 1

fi

# Holds path to script directory
script_path=$(dirname "$script_path")

# Creates working directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

# Creates temp directory in working directory if none exist
temp="$WD/temp"
[ -d "$temp" ] || mkdir -m 775 "$temp"

# Creates data directory in working directory if none exist
data="$WD/03_sfs"
[ -d "$data" ] || mkdir -m 775 "$data"

# Creates a species specific sample directory within the data directory
samplename="$(basename "$(dirname "$S")")"
speciesname="$(basename "$(dirname "$(dirname "$S")")")"
sampledir="$data/$speciesname/$samplename"
[ -d "$sampledir" ] || mkdir -p "$sampledir"
chmod 775 "$data/$speciesname" "$sampledir"

# Creates directory for .out files
out="$sampledir"/out
[ -d "$out" ] || mkdir -p "$out"
chmod 775 "$out"

# Creates file with all QNAME from .bam file
qname="$temp/$(basename "$(dirname "$S")")_qname.txt"
if [ ! -e "$qname" ]; then

    samtools view -H "$(dirname "$S")"/"$(basename "$(dirname "$S")")"_filtered.bam \
    | grep '@SQ' \
    | awk '{print $2}' \
    > "$temp"/"$(basename "$(dirname "$S")")"_headertemp.txt

    while read -r line; do
        echo "${line:3}"
    done < "$temp"/"$(basename "$(dirname "$S")")"_headertemp.txt > "$qname"

    # Removes temp header file
    rm -f "$temp"/"$(basename "$(dirname "$S")")"_headertemp.txt

fi

# Number of lines, eg. number of unique IDs
lines=$(wc -l <"$qname")

# Function to Calculates time adjustment based on filesize and number of IDs comparative to .mpileup of Ocin_NYS-F
filesize=$(wc -c <"$S")
sizeratio=$(awk -v filesize="$filesize" 'BEGIN { print ( filesize / 234060585564) }')
lineratio=$(awk -v lines="$lines" 'BEGIN { print ( 9402 / lines ) }')
adjustment=$(awk -v sizeratio="$sizeratio" -v lineratio="$lineratio" 'BEGIN { print ( sizeratio * lineratio ) }')

# Checks whether to force owerwrite
if [ "$force_overwrite" == "Y" ]; then
    runlist=("0" "0" "0" "0" "0")
else
    filechecker
fi

job1 "$lines" "$qname" "$S" "$A" "$namebase" "$temp" "$sampledir" "$keep_temp"

job2 "$lines" "$qname" "$S" "$A" "$namebase" "$temp" "$sampledir" "$keep_temp"

job3 "$lines" "$qname" "$S" "$A" "$namebase" "$temp" "$sampledir" "$keep_temp"

job4 "$lines" "$qname" "$S" "$A" "$namebase" "$temp" "$sampledir" "$keep_temp"

job5 "$lines" "$qname" "$S" "$A" "$namebase" "$temp" "$sampledir" "$keep_temp"

job6 "$lines" "$qname" "$S" "$A" "$namebase" "$temp" "$sampledir" "$keep_temp"

# Clean up of empty stdout files
sbatch \
	--output=/dev/null \
	--error=/dev/null \
	--dependency=afterany:"$jid6":"$jid7":"$jid8" \
	"$script_path"/03_cleanup.sh "$lines" "$qname" "$S" "$A" "$namebase" "$temp" "$sampledir" "$keep_temp"

exit 0