#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 04_assembly_init.sh kmer_analysis -s <FILE> -k <STRING> [OPTIONS]

This script utilizes jellyfish and GenomeScope to conduct a kmer analysis
on a genome file in fastx format.

OPTIONS:
    -s  FILE        Target file containing raw genome sequence.
    -k  STRING      Kmer interval in the format: '[start]-[end]:[step]',
                    '[start]' being an integer signifying the first kmer length,
                    '[end]' being an integer signifying the last kmer length,
                    and '[step]', being an integer signifying the step-size.
                    If a single kmer length is desired write '[start]-[end]',
                    where '[start]' and '[end]' are identical.
    -c              Use cannonical kmer count. Need to be set when running on
                    sequencing reads. If running on an actual genome or finished
                    sequence do NOT set.
    -m  INTEGER     Cut-off for excluding high frequency kmers from analysis.
                    Default is '1000'.
    -d  DIRECTORY   Desired working directory. Defaults to current directory.

    -h              Show this message

Original:
    https://github.com/gmarcais/Jellyfish
    https://github.com/schatzlab/genomescope
    
EOF
}

    # -l  INTEGER     Read length. Default is '100'. (Only used for GenomeScope1)
    # -2              Switches from use of GenomeScope1 to GenomeScope2.

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

# Default cut-off value for max kmer
kmer_max="1000"

# Default read length value
read_length="100"

# Used for switching environments
currentuser="jepe"

# Default value for canonical option
canonical=

# ----------------- Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 's:k:cm:d:2hu:' OPTION; do
    case "$OPTION" in
		s)
	    	if [ -f "$OPTARG" ]; then
                targetfile="$(readlink -f "$OPTARG")"
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
				exit 1 
            fi
	    	;;
        k)
            kmer_interval="$OPTARG"
            ;;
        c)
            canonical="1"
            ;;
        m)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                kmer_max="$OPTARG"
            else
                echo -e "\n-m can only be assinged an integer value\n"
                exit 1
            fi
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
    WD="$WD/kmer_analysis"
    [ -d "$WD" ] || mkdir -m 775 "$WD"
else
    WD="$(dirname "$targetfile")/kmer_analysis"
    [ -d "$WD" ] || mkdir -m 775 "$WD"
fi

out="$WD/out"
[ -d "$out" ] || mkdir -m 775 "$out"

jid1=$(sbatch \
    --parsable \
    --array="$kmer_interval" \
    --time=240 \
    --mem=480G \
    --cpus-per-task=32 \
    --output="$out"/kmer_count-%a-%j.out \
    "$script_path"/modules/04_kmer_01_count.sh "$targetfile" "$WD" "$currentuser" "$canonical")

jid2=$(sbatch \
    --parsable \
    --array="$kmer_interval" \
    --time=240 \
    --mem=480G \
    --cpus-per-task=32 \
    --dependency=aftercorr:"$jid1" \
    --output="$out"/kmer_histogram-%a-%j.out \
    "$script_path"/modules/04_kmer_02_histogram.sh "$targetfile" "$WD" "$currentuser")

# shellcheck disable=2034
# jid3=$(sbatch \
#     --parsable \
#     --array="$kmer_interval"\
#     --time=30 \
#     --mem=20G \
#     --cpus-per-task=2 \
#     --dependency=aftercorr:"$jid2" \
#     --output="$out"/kmer_genomescope-%a-%j.out \
#     "$script_path"/modules/04_kmer_03_genomescope.sh "$targetfile" "$WD" "$script_path" "$kmer_max" "$read_length" "$currentuser")

# shellcheck disable=2034
jid3=$(sbatch \
    --parsable \
    --array="$kmer_interval"\
    --time=240 \
    --mem=20G \
    --cpus-per-task=2 \
    --dependency=aftercorr:"$jid2" \
    --output="$out"/kmer_genomescope2-%a-%j.out \
    "$script_path"/modules/04_kmer_03_genomescope2.sh "$targetfile" "$WD" "$script_path" "$kmer_max" "$read_length" "$currentuser")


exit 0