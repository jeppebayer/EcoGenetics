#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 04_assembly_init.sh hifiadapterfilt -s <FILE>

This script finds and removes remaining adapters specifically from
PacBio HiFi reads.

OPTIONS:
    -s  FILE        Target file containing raw PacBio HiFi reads

    -h              Show this message

OUTPUT:
    {input_name}.contaminant.blastout   (Output from BLAST search)
    {input_name}.blocklist              (Headers of PacBio adapter
                                        contaminated reads to be removed)
    {input_name}.filt.fastq.gz          (Fastq reads free of PacBio
                                        adapter sequence ready for 
                                        assembly)
    {input_name}.stats                  (File with simple math on number
                                        of reads removed, etc)

Original:
    https://github.com/sheinasim/HiFiAdapterFilt
    
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
                targetfile="$(readlink -f "$OPTARG")"
                targetdir="$(dirname "$targetfile")"
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

# Creates directories
filtdir="$targetdir"/HiFiAdapterFilt
[ -d "$filtdir" ] || mkdir -m 775 "$filtdir"
out="$filtdir"/out
[ -d "$out" ] || mkdir -m 775 "$out"

# Checks file extension and creates correct file prefix
if [[ "$targetfile" == *.bam ]]; then
    targetprefix=${targetfile%.*}
    ext="bam"
elif [[ "$targetfile" == *.fastq ]]; then
    targetprefix=${targetfile%.*}
    ext="fastq"
elif [[ "$targetfile" == *.fastq.gz ]]; then
    targetprefix=${targetfile%.*.*}
    ext="fastq.gz"
elif [[ "$targetfile" == *.fq ]]; then
    targetprefix=${targetfile%.*}
    ext="fq"
elif [[ "$targetfile" == *.fq.gz ]]; then
    targetprefix=${targetfile%.*.*}
    ext="fq.gz"
fi

# Checks if other files in the folder share the same prefix
filecount=0
if [ -f "$targetprefix".bam ]; then
    ((filecount++))
fi
if [ -f "$targetprefix".fastq ]; then
    ((filecount++))
fi
if [ -f "$targetprefix".fastq.gz ]; then
    ((filecount++))
fi
if [ -f "$targetprefix".fq ]; then
    ((filecount++))
fi
if [ -f "$targetprefix".fq.gz ]; then
    ((filecount++))
fi

# If more than 1 file has the target prefix creates a copy with a unique prefix
if [ $((filecount)) -gt 1 ]; then
    cp "$targetfile" "$targetdir"/temp_name."$ext"
    originalprefix="$targetprefix"
    targetprefix="temp_name"
fi

# shellcheck disable=2034
jid1=$(sbatch \
	--parsable \
    --chdir="$targetdir" \
	--time=120 \
	--mem-per-cpu=10G \
	--cpus-per-task=10 \
	--output="$out"/04_hifiadapterfilt-%j.out \
	"$script_path"/modules/04_hifiadapterfilt_01_hifiadapterfilt.sh "$targetprefix" "$filtdir" "$script_path" "$currentuser" "$originalprefix")

exit 0