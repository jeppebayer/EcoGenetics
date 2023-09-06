#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 04_assembly_init.sh blobtools -s <FILE> -a <FILE> [OPTIONS]

This script runs blobtools, which utilizes coverage and sequence blasting
(NCBI BLAST+ and DIAMOND), to remove possible contaminants from genome assembly.

OPTIONS:
    -s  FILE        File containing raw HiFI reads.
    -a  FILE        File containg genome assembly corresponding to selected
                    raw HiFi reads.
    -d  DIRECTORY   Desired working directory. Defaults to current directory.
    -b  INTERGER    Start from BLAST. Argument can be '1' for ncbimegablast, '2'
                    for diamondblast, or '12' | '21' for both. (This option is
                    generally only to be used for troubleshooting)

    -h              Show this message

Original:
    https://github.com/DRL/blobtools
    https://github.com/bbuchfink/diamond
    
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

# Switch to start from BLAST
startfromblast="N"

# ----------------- Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 's:a:d:b:hu:' OPTION; do
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
        b)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                startfromblast="Y"
                whichblast="$OPTARG"
            else
                echo -e "\n$OPTARG is an invalid argument to -b\n"
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


WD="$WD/04_assembly"
[ -d "$WD" ] || mkdir -m 775 "$WD"

if [ "$(basename "$(dirname "$rawreads")")" == "HiFiAdapterFilt" ]; then
    species_name="$(basename "$(dirname "$(dirname "$rawreads")")")"
else
    species_name="$(basename "$(dirname "$rawreads")")"
fi
WD="$WD"/"$species_name"
[ -d "$WD" ] || mkdir -m 775 "$WD"  
WD="$WD/blobtools"
[ -d "$WD" ] || mkdir -m 775 "$WD"
out="$WD/out"
[ -d "$out" ] || mkdir -m 775 "$out"

# Creates species abbreviation from species name
genus=${species_name%_*}; genus=${genus::3}; genus=${genus^}
species=${species_name#*_}; species=${species::3}; species=${species^}
speciesabbr="$genus""$species"

if [ "$startfromblast" == "N" ]; then
    jid1=$(sbatch \
        --parsable \
        --time=180 \
        --mem-per-cpu=15G \
        --cpus-per-task=16 \
        --output="$out"/align-%j.out \
        "$script_path"/modules/04_blobtools_01_align.sh "$rawreads" "$asmfile" "$WD" "$temp" "$speciesabbr" "$currentuser")

    jid2=$(sbatch \
        --parsable \
        --time=180 \
        --mem-per-cpu=10G \
        --cpus-per-task=2 \
        --dependency=afterany:"$jid1" \
        --output="$out"/cov-%j.out \
        "$script_path"/modules/04_blobtools_02_cov.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

    jid3_1=$(sbatch \
        --parsable \
        --time=8640 \
        --mem-per-cpu=15G \
        --cpus-per-task=16 \
        --output="$out"/blastn-%j.out \
        "$script_path"/modules/04_blobtools_03_blastn.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

    jid3_2=$(sbatch \
        --parsable \
        --time=8640 \
        --mem-per-cpu=15G \
        --cpus-per-task=16 \
        --output="$out"/blastd-%j.out \
        "$script_path"/modules/04_blobtools_03_blastd.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

    jid4=$(sbatch \
        --parsable \
        --time=60 \
        --mem-per-cpu=10G \
        --cpus-per-task=2 \
        --dependency=afterany:"$jid2":"$jid3_1":"$jid3_2" \
        --output="$out"/create-%j.out \
        "$script_path"/modules/04_blobtools_04_create.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

elif [ "$startfromblast" == "Y" ]; then
    if [ "$whichblast" == "1" ]; then
        jid3_1=$(sbatch \
            --parsable \
            --time=8640 \
            --mem-per-cpu=15G \
            --cpus-per-task=16 \
            --output="$out"/blastn-%j.out \
            "$script_path"/modules/04_blobtools_03_blastn.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")
        
        jid4=$(sbatch \
            --parsable \
            --time=60 \
            --mem-per-cpu=10G \
            --cpus-per-task=2 \
            --dependency=afterany:"$jid3_1" \
            --output="$out"/create-%j.out \
            "$script_path"/modules/04_blobtools_04_create.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

    elif [ "$whichblast" == "2" ]; then
        jid3_2=$(sbatch \
            --parsable \
            --time=8640 \
            --mem-per-cpu=15G \
            --cpus-per-task=16 \
            --output="$out"/blastd-%j.out \
            "$script_path"/modules/04_blobtools_03_blastd.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")
        
        jid4=$(sbatch \
            --parsable \
            --time=60 \
            --mem-per-cpu=10G \
            --cpus-per-task=2 \
            --dependency=afterany:"$jid3_2" \
            --output="$out"/create-%j.out \
            "$script_path"/modules/04_blobtools_04_create.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

    else
        jid3_1=$(sbatch \
            --parsable \
            --time=8640 \
            --mem-per-cpu=15G \
            --cpus-per-task=16 \
            --output="$out"/blastn-%j.out \
            "$script_path"/modules/04_blobtools_03_blastn.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

        jid3_2=$(sbatch \
            --parsable \
            --time=8640 \
            --mem-per-cpu=15G \
            --cpus-per-task=16 \
            --output="$out"/blastd-%j.out \
            "$script_path"/modules/04_blobtools_03_blastd.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

        jid4=$(sbatch \
            --parsable \
            --time=60 \
            --mem-per-cpu=10G \
            --cpus-per-task=2 \
            --dependency=afterany:"$jid3_1":"$jid3_2" \
            --output="$out"/create-%j.out \
            "$script_path"/modules/04_blobtools_04_create.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")
    fi
fi

# jid4=$(sbatch \
#     --parsable \
#     --time=600 \
#     --mem-per-cpu=10G \
#     --cpus-per-task=2 \
#     --dependency=afterany:"$jid3_1":"$jid3_2" \
#     --output="$out"/create-%j.out \
#     "$script_path"/modules/04_blobtools_04_create.sh "$asmfile" "$WD" "$speciesabbr" "$currentuser")

# # shellcheck disable=2034
# jid5=$(sbatch \
#     --parsable \
#     --time=180 \
#     --mem-per-cpu=10G \
#     --cpus-per-task=2 \
#     --dependency=afterany:"$jid4" \
#     --output="$out"/plot-%j.out \
#     "$script_path"/modules/04_blobtools_05_plot.sh "$WD" "$speciesabbr" "$currentuser")

# # shellcheck disable=2034
# jid6=$(sbatch \
#     --parsable \
#     --time=120 \
#     --mem-per-cpu=10G \
#     --cpus-per-task=2 \
#     --dependency=afterany:"$jid4" \
#     --output="$out"/view-%j.out \
#     "$script_path"/modules/04_blobtools_06_view.sh "$WD" "$speciesabbr" "$currentuser")

exit 0