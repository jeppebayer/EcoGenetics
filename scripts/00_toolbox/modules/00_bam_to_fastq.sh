#!/bin/bash

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate data_prep
fi

run_jobs=false

# ------------------ Flag Processing --------------------------------------

while getopts 's:d:rb:' OPTION; do
    case "$OPTION" in
        s)
            if [ -d "$OPTARG" ]; then
                species_dir=$(readlink -f "$OPTARG")
            else
                echo -e "\n-s has to be supplied a directory.\n"
            fi
            ;;
        d)
            if [ -d "$OPTARG" ]; then
                target_dir=$(readlink -f "$OPTARG")
            else
                echo -e "\n-d has to be an existing directory.\n"
                exit 1
            fi
            ;;
        r)
            run_jobs=true
            ;;
        b)
            unmappedbam="$OPTARG"
            ;;
        ?)
            exit 2
            ;;
    esac
done

# Removes "/" from the end of directory path
if [ "${species_dir:-1}" == "/" ]; then
    length=${#species_dir}
    species_dir=${species_dir::length - 1}
fi

# Removes "/" from the end of directory path
if [ "${target_dir:-1}" == "/" ]; then
    length=${#target_dir}
    target_dir=${target_dir::length - 1}
fi

# ------------------ Main -------------------------------------------------

if [ "$run_jobs" == true ]; then
    target_dir="$target_dir"/"$(basename "$species_dir")"
    # Creates specified directory if it doesn't exist
    [ -d "$target_dir" ] || mkdir -m 775 "$target_dir"
    sample_name=$(basename "$unmappedbam")
    sample_name=${sample_name%.*}
    samtools fastq -@ 6 -0 /dev/null "$unmappedbam" > "$target_dir"/"$sample_name".fq.gz

elif [ "$run_jobs" == false ]; then
    counter=0
    jidarray=()
    for unmappedbam in "$species_dir"/*/*_unmapped.bam; do
        sample_name=$(basename "$unmappedbam")
        sample_name=${sample_name%.*}
        if [ ! -f "$target_dir"/"$(basename "$species_dir")"/"$sample_name".fq.gz ]; then
            jid=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir=/faststorage/project/EcoGenetics \
                --time=30 \
                --mem=36G \
                --cpus-per-task=6 \
                --output=/dev/null \
                --error=/dev/null \
                "$(readlink -f "$0")" -s "$species_dir" -d "$target_dir" -b "$unmappedbam" -r)
            ((counter++))
            jidarray+=("$jid")
        else
            echo -e "\n$sample_name.bam is skipped. FASTQ version already exists"
        fi
    done
    echo -e "\nA total of $counter jobs have been sent to the queue\nThe associated JobIDs are as follows: ${jidarray[*]}\nFASTQ files can be found in $target_dir/$(basename "$species_dir")/\n"
fi

exit 0