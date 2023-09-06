#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 04_assembly_init.sh purge_dups -s <FILE> -a <FILE> [OPTIONS]

This script runs purge_dups, in combination with minimap2, attempting
to remove haplotigs and overlaps from a genome assembly based on read depth.

OPTIONS:
    -s  FILE        File containing raw HiFI reads.
    -a  FILE        File containg genome assembly corresponding to selected
                    raw HiFi reads.
    -d  DIRECTORY   Desired working directory. Defaults to current directory.

    -h              Show this message

Original:
    https://github.com/dfguan/purge_dups
    https://github.com/lh3/minimap2
    
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

# Used for switching environments
currentuser="jepe"

nonstandard=false

# ----------------- Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 's:a:d:hu:x' OPTION; do
    case "$OPTION" in
		s)
	    	if [ -f "$OPTARG" ]; then
                rawreads="$(readlink -f "$OPTARG")"
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
				exit 1 
            fi
	    	;;
        a)
	    	if [ -f "$OPTARG" ]; then
                asmfile="$(readlink -f "$OPTARG")"
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
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
        x)
            nonstandard=true
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

if [ "$nonstandard" = false ]; then
    WD="$WD/04_assembly"
    [ -d "$WD" ] || mkdir -m 775 "$WD"

    if [ "$(basename "$(dirname "$rawreads")")" == "HiFiAdapterFilt" ]; then
        species_name="$(basename "$(dirname "$(dirname "$rawreads")")")"
    else
        species_name="$(basename "$(dirname "$rawreads")")"
    fi
    genus=${species_name%_*}
    genus=${genus::3}
    genus=${genus^}
    species=${species_name#*_}
    species=${species::3}
    species=${species^}
    speciesabbr="$genus""$species"

    WD="$WD"/"$species_name"
    [ -d "$WD" ] || mkdir -m 775 "$WD"
else
    if [ "$(basename "$(dirname "$(dirname "$asmfile")")")" == "purge_dups" ]; then
        WD="$(dirname "$(dirname "$asmfile")")"
    else
        WD="$(dirname "$asmfile")"
        WD="$WD/purge_dups"
        [ -d "$WD" ] || mkdir -m 775 "$WD"
    fi
    name="$(basename "$rawreads")"
    oldIFS=$IFS
    IFS='.'
    read -ra name <<< "$name"
    IFS=$oldIFS
    speciesabbr=${name[0]}
fi

count=$(find "$WD" -maxdepth 1 -mindepth 1 -type d | wc -l)
count=$((count + 1))
purge_num="0$count"
WD="$WD/$purge_num"
[ -d "$WD" ] || mkdir -m 775 "$WD"
out="$WD/out"
[ -d "$out" ] || mkdir -m 775 "$out"

jid1=$(sbatch \
    --parsable \
    --time=60 \
    --mem-per-cpu=10G \
    --cpus-per-task=15 \
    --output="$out"/minimap2_1-%j.out \
    "$script_path"/modules/04_purge_dups_01_minimap2.sh "$rawreads" "$asmfile" "$WD" "$speciesabbr" "$currentuser")

jid2=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid1" \
    --output="$out"/pbcstat-%j.out \
    "$script_path"/modules/04_purge_dups_02_pbcstat.sh "$WD" "$speciesabbr" "$currentuser")

jid3=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid2" \
    --output="$out"/calcuts-%j.out \
    "$script_path"/modules/04_purge_dups_03_calcuts.sh "$WD" "$currentuser")

jid4=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid3" \
    --output="$out"/split-%j.out \
    "$script_path"/modules/04_purge_dups_04_split.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

jid5=$(sbatch \
    --parsable \
    --time=60 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid4" \
    --output="$out"/minimap2_2-%j.out \
    "$script_path"/modules/04_purge_dups_05_minimap2.sh "$WD" "$speciesabbr" "$currentuser")

jid6=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=15 \
    --dependency=afterany:"$jid5" \
    --output="$out"/purge-%j.out \
    "$script_path"/modules/04_purge_dups_06_purge.sh "$WD" "$speciesabbr" "$currentuser")

# shellcheck disable=2034
jid7=$(sbatch \
    --parsable \
    --chdir="$WD" \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid6" \
    --output="$out"/get_seqs-%j.out \
    "$script_path"/modules/04_purge_dups_07_get_seqs.sh "$asmfile" "$WD" "$speciesabbr" "$script_path" "$currentuser")

exit 0
