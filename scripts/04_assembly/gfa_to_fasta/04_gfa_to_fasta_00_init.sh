#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 04_assembly_init.sh gfa_to_fasta -s <FILE> [OPTIONS]

Small script to convert .gfa file to simple fasta format.

OPTIONS:
    -s  FILE        File in .gfa format.

    -h              Show this message
    
EOF
}

# ----------------- Configuration ----------------------------------------

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

# Used for switching environments
currentuser="jepe"

# ----------------- Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 's:hu:' OPTION; do
    case "$OPTION" in
		s)
	    	if [ -f "$OPTARG" ]; then
                if [[ "$OPTARG" == *.gfa ]]; then
                    targetfile="$(readlink -f "$OPTARG")"
                else
                    echo -e "\n$OPTARG is not in .gfa format\n"
				    exit 1
                fi
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
				exit 1 
            fi
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

# ----------------- Script Queue -----------------------------------------

# shellcheck disable=2034
jid1=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --output=/dev/null \
	--error=/dev/null \
    "$script_path"/modules/04_gfa_to_fasta_01_gfa_to_fasta.sh "$targetfile" "$currentuser")

exit 0