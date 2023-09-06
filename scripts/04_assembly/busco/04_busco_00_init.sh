#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 04_assembly_init.sh busco -s <FILE> [OPTIONS]

This script runs a BUSCO analyses on a genome assembly from HiFi reads.

OPTIONS:
    -s  FILE        Target file containing genome assembly.
    -d  DIRECTORY   Desired working directory. Defaults to current directory.
    -g  STRING      If 'list' will show available datasets which can be passed
                    for the busco analysis. arthropoda_odb10 [default].
    -t  INTEGER     Extra time to add in hours. Use if time calculation is
                    expected to be low.

    -h              Show this message

Original:
    https://busco.ezlab.org/
    
EOF
}

# ----------------- Configuration ----------------------------------------

# Working directory
WD="$(readlink -f "$PWD")"

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    # shellcheck disable=2153
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

script_path="$(dirname "$script_path")"

# Standard dataset & dataset list
dataset="arthropoda_odb10"
datasetlist="$(ls /faststorage/project/EcoGenetics/BACKUP/database/busco/busco_downloads/lineages/)"

# Standard additional time
addtime=0

# Used for switching environments
currentuser="jepe"

# ----------------- Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 's:d:g:t:hu:' OPTION; do
    case "$OPTION" in
		s)
	    	if [ -f "$OPTARG" ]; then
                targetfile="$(readlink -f "$OPTARG")"
                targetdir="$(dirname "$targetfile")"
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
				exit 1 
            fi
	    	;;
        d)
            WD="$(readlink -f "$OPTARG")"  
            ;;
        g)
            if [ "$OPTARG" == "list" ]; then
                echo -e "\n$datasetlist\n"
                exit 0
            else
                dataset="$OPTARG"
            fi
            ;;
        t)
            addtime="$OPTARG"
            addtime=$((addtime * 60))
            ;;
        h)
            usage
            exit 0
            ;;
        u)
            currentuser="$OPTARG"
            ;;
        :)
            echo "Must supply an argument to -$OPTARG"
            exit 1
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

# ----------------- Functions --------------------------------------------

# Function used to adjust time for jobs
timer()
{
    jobtime=$(awk -v filesize="$filesize" -v addtime=$addtime 'BEGIN {print int( ((306 * log(filesize)) - 5680) * 1.5 + addtime)}')
    echo "$jobtime"
}

# ----------------- Script Queue -----------------------------------------

# Creates working directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

temp="$WD"/temp
[ -d "$temp" ] || mkdir -m 775 "$temp"

if [[ "$targetfile" == *BACKUP* ]]; then
    WD="$WD"/04_assembly
    [ -d "$WD" ] || mkdir -m 775 "$WD"
    WD="$WD"/"$(basename "$(dirname "$targetfile")")"
    [ -d "$WD" ] || mkdir -m 775 "$WD"
else
    WD="$targetdir"
    [ -d "$WD" ] || mkdir -m 775 "$WD"
fi

WD="$WD"/busco
[ -d "$WD" ] || mkdir -m 775 "$WD"

out="$WD/out"
[ -d "$out" ] || mkdir -m 775 "$out"

filesize=$(wc -c < "$targetfile")

# shellcheck disable=2034
jid1=$(sbatch \
    --parsable \
    --chdir="$WD" \
    --time="$(timer)" \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --output="$out"/busco-%j.out \
    "$script_path"/modules/04_busco_01_busco.sh "$targetfile" "$WD" "$currentuser" "$dataset")

exit 0