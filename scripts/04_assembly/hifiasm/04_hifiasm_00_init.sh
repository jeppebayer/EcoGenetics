#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 04_assembly_init.sh hifiasm -s <FILE> [OPTIONS]

This script utilizes hifiasm to create a genome assembly from HiFi reads.

OPTIONS:
    -s  FILE        Target file containing raw PacBio HiFi reads.
    -t  INTEGER     Number of bases to trim off each end of all sequences.
    -l  FLOAT       Similarity threshold for duplicate haplotigs. 
                    Default is '0.55'. If sample has high heterozygosity,
                    then lower the value, start with '0.1'.
    -d  DIRECTORY   Desired working directory. Defaults to current directory.

    -h              Show this message

Original:
    https://github.com/chhylp123/hifiasm
    
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

# Default trim value
trim=

# Default similarity threshold value
sim_thres=0.55

# Used for switching environments
currentuser="jepe"

# ----------------- Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 's:t:l:d:hu:' OPTION; do
    case "$OPTION" in
		s)
	    	if [ -f "$OPTARG" ]; then
                targetfile="$(readlink -f "$OPTARG")"
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
				exit 1 
            fi
	    	;;
        t)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                trim="$OPTARG"
            else
                echo -e "\n-t can only be assinged an integer value"
            fi
            ;;
        l)
            sim_thres="$OPTARG"
            ;;
        d)
            WD="$(readlink -f "$OPTARG")"  
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

# ----------------- Script Queue -----------------------------------------

# Creates working directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

temp="$WD/temp"
[ -d "$temp" ] || mkdir -m 775 "$temp"

if [[ "$targetfile" == *BACKUP* ]]; then
    WD="$WD/04_assembly"
    [ -d "$WD" ] || mkdir -m 775 "$WD"
    if [ "$(basename "$(dirname "$targetfile")")" == "HiFiAdapterFilt" ]; then
        species_name="$(basename "$(dirname "$(dirname "$targetfile")")")"
        WD="$WD/$species_name"
    else
        species_name="$(basename "$(dirname "$targetfile")")"
        WD="$WD/$species_name"
    fi
    [ -d "$WD" ] || mkdir -m 775 "$WD"
    genus=${species_name%_*}
    genus=${genus::3}
    genus=${genus^}
    species=${species_name#*_}
    species=${species::3}
    species=${species^}
    speciesabbr="$genus""$species"   
else
    WD="$(dirname "$targetfile")/hifiasm"
    [ -d "$WD" ] || mkdir -m 775 "$WD"
    species_name="$(basename "$targetfile")"
    speciesabbr=${species_name:0:6}
fi

out="$WD/out"
[ -d "$out" ] || mkdir -m 775 "$out"

# shellcheck disable=2034
jid1=$(sbatch \
    --parsable \
    --time=600 \
    --mem-per-cpu=15G \
    --cpus-per-task=32 \
    --output="$out"/hifiasm-%j.out \
    "$script_path"/modules/04_hifiasm_01_hifiasm.sh "$targetfile" "$WD" "$trim" "$sim_thres" "$speciesabbr" "$currentuser")

exit 0