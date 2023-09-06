#!/bin/bash

# ----------------- Configuration ----------------------------------------

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate gtdb
fi

shopt -s extglob

send_to_queue=true

# ----------------- Usage ------------------------------------------------

usage(){
cat << EOF

Usage: 00_gtdb_accession2taxid.sh -d <FILE> -n <FILE> -a <FILE> [OPTIONS]

Uses the names.dmp file associated with NCBI database to create an accession2taxid
file which links GTDB accession numbers to NCBI taxIDs.

    -d  FILE            GTDB file in fasta format containing accession number,
                        species name and sequence.
    -n  FILE            NCBI names.dmp file.
    -a  FILE            NCBI accession2taxid file.

OPTIONS:
    -h                  Shows this message.

EOF
}

# ------------------ Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 'd:n:a:qh' OPTION; do
    case "$OPTION" in
        d)
            if [[ "$OPTARG"  == *.@(fna|fa|fasta) ]]; then
                gtdb_seq=$(readlink -f "$OPTARG")
                target_dir=$(dirname "$gtdb_seq")
            else
                echo -e "\nFile supplied to -d has to be in one of the following formats: .fna, .fa, .fasta\n"
                exit 1
            fi
            ;;
        n)
            if [[ "$OPTARG" == *.dmp ]]; then
                name_dump=$(readlink -f "$OPTARG")
            else
                echo -e "\nFile supplied to -n has to be a .dmp\n"
                exit 1
            fi
            ;;
        a)
            if [[ "$OPTARG" == *.accession2taxid* ]]; then
                ncbi_accession=$(readlink -f "$OPTARG")
            else
                echo -e "\nFile supplied to -a must be named in the format: '*.accession2taxid*'"
            fi
            ;;
        q)
            send_to_queue=false
            ;;
        h)
            usage
            exit 0
            ;;
        ?)
            usage
            exit 2
            ;;
    esac
done

# Tests that GTDB file has been supplied
if [[ -z $gtdb_seq ]]; then
    echo -e "\n-d must be assigned\n"
    exit 3
fi

# Tests that names.dmp has been supplied
if [[ -z $name_dump ]]; then
    echo -e "\n-n must be assigned\n"
    exit 3
fi

# Tests that NCBI accession2taxid file has been supplied
if [[ -z $ncbi_accession ]]; then
    echo -e "\n-a must be assigned\n"
    exit 3
fi

# ----------------- Functions --------------------------------------------

# ------------------ Main -------------------------------------------------

if [ "$send_to_queue" == true ]; then
    jid=$(sbatch \
        --parsable \
        --account=EcoGenetics \
        --chdir="$target_dir" \
        --time=7-00:00:00 \
        --mem=48G \
        --cpus-per-task=8 \
        --output="$target_dir"/gtdb_accession2taxid-%j.out \
        "$(readlink -f "$0")" -d "$gtdb_seq" -n "$name_dump" -a "$ncbi_accession" -q)

elif [ "$send_to_queue" == false ]; then
    echo -e "JobID:\t$SLURM_JOB_ID"
    start_date="$(date +"%Y-%m-%d")"
    start_time=$(date +"%T")
    timerequested=$(squeue -j "$SLURM_JOB_ID" -h --Format TimeLimit)
    echo -e "Start date:\t$start_date\nStart time:\t$start_time\nTime allocated to job:\t$timerequested"
    echo -e "\nGTDB sequence file used:\t$gtdb_seq\nNCBI names.dmp file used:\t$name_dump\nNCBI accession2taxid file used:\t$ncbi_accession"

    # Gets name of GTDB file without extension
    gtdb_seq_name=$(basename "$gtdb_seq"); gtdb_seq_name=${gtdb_seq_name%.*}
    
    # ============================== Part 1 ==============================
    # Creates file containing only headers from GTDB file
    headerfile="$target_dir"/"$gtdb_seq_name"_headers.txt
    echo -e "\nCreating file containing all headers from GTDB sequence file..."
    SECONDS=0
    bioawk -c fastx '{printf ("%s\t%s\n", $name, $comment)}' "$gtdb_seq" > "$headerfile"
    duration=$SECONDS
    echo -e "Process took:\t$((duration / 3600)) hours $((duration / 60 % 60)) minutes $((duration % 60)) seconds"
    
    # Gets number of lines in headerfile
    nentries=$(wc -l < "$headerfile")
    echo -e "Number of entires in GTDB sequence file:\t$nentries"

    # ============================== Part 2 ==============================
    ncbi_accession_dir="$target_dir"/ncbi_accession2taxid_parts
    # Creates specified directory if it doesn't exist
    [ -d "$ncbi_accession_dir" ] || mkdir -m 775 "$ncbi_accession_dir"
    echo -e "\nDirectory containing alphabetically split parts of NCBI accession2taxid file is located:\t$ncbi_accession_dir"
    
    # Creates separate files for each accession ID starting letter in NCBI accession2taxID file
    echo -e "\nSplitting the NCBI accession2taxid file..."
    SECONDS=0
    alphabet=(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z)
    # Checks whether or not the NCBI accession2taxid file is compressed
    if [[ "$ncbi_accession" == *.gz ]]; then
        ncbi_accession_name=$(basename "${ncbi_accession%.*}")
        for i in "${alphabet[@]}"; do
            head -n 1 <(zcat "$ncbi_accession") > "$ncbi_accession_dir"/"$i"_"$ncbi_accession_name"
            awk -v letter="$i" 'index($1, letter)==1' <(zcat "$ncbi_accession") >> "$ncbi_accession_dir"/"$i"_"$ncbi_accession_name"
        done
    else
        ncbi_accession_name=$(basename "$ncbi_accession")
        for i in "${alphabet[@]}"; do
            head -n 1 "$ncbi_accession" > "$ncbi_accession_dir"/"$i"_"$ncbi_accession_name"
            awk -v letter="$i" 'index($1, letter)==1' "$ncbi_accession" >> "$ncbi_accession_dir"/"$i"_"$ncbi_accession_name"
        done
    fi
    duration=$SECONDS
    echo -e "Process took:\t$((duration / 3600)) hours $((duration / 60 % 60)) minutes $((duration % 60)) seconds"

    # ============================== Part 3 ==============================
    # Finds GTDB accession IDs already in NCBI data by searching only in NCBI accession2taxid file parts relevant
    echo -e "\nChecking for duplicate accession IDs..."
    SECONDS=0
    duplicate_array=()
    for ((i=1; i<=nentries; i++)); do
        starting_character="$(awk -v linenum="$i" 'NR == linenum {print substr($1,0,1)}' "$headerfile")"
        if rg -j 8 -Fq <(awk -v linenum="$i" 'NR == linenum {print $1}' "$headerfile") <(awk '{print $2}' "$ncbi_accession_dir"/"$starting_character"_"$ncbi_accession_name"); then
            duplicate_array+=("$i")
            # echo -e "Duplicate found at line:\t$i"
        fi
        tcounter=1
        thousands=$((i / 1000))
        if [ $thousands -eq $tcounter ]; then
            echo -e -n "."
            ((tcounter++))
        fi
    done
    echo 
    for entry in "${duplicate_array[@]}"; do
        echo -e "Duplicate found at line:\t$entry"
    done
    # Removes lines from headerfile which contains accession IDs already in NCBI file.
    # Starts from the bottom of the file as not change relevant line numbers for subsequent entries
    for ((i=${#duplicate_array[@]} - 1; i>=0; i--)); do
        line_number=${duplicate_array[$i]}
        # echo "$line_number" >> "$target_dir"/duplicate_line_num.txt
        sed -i "{$line_number}d" "$headerfile"
    done
    duration=$SECONDS
    echo -e "Process took:\t$((duration / 3600)) hours $((duration / 60 % 60)) minutes $((duration % 60)) seconds"
    echo -e "Number of duplicates found:\t${#duplicate_array[@]}\n"

    # Gets new number of lines in headerfile
    nentries=$(wc -l < "$headerfile")

    # ============================== Part 4 ==============================
    # Creates gtdb_accession2taxid file and inserts header
    gtdb_accession2taxid_file="$target_dir"/gtdb.accession2taxid
    echo -e "accession\taccession.version\taxid\tgi" > "$gtdb_accession2taxid_file"

    echo -e "\nExtracting taxIDs matching species names for accession IDs..."
    SECONDS=0
    entrynum=0
    # Loops through each entry in headerfile
    for ((i=1; i<=nentries; i++)); do
        # Gets comments related to entry, tries to isolate species name
        original_species_name="$(awk -v linenum="$i" 'NR == linenum {print $2}' "$headerfile")"
        # Removes first instance of ',' from the right and anything that comes after
        species_name=${original_species_name%,*}
        # Removes first instance of 'isolate' from the left and anything that comes after
        species_name=${species_name%% isolate*}
        # Removes first instance of ':' from the left side and anything that comes before
        species_name=${species_name#*: }
        # Removes first instance of 'uncultured' from the left side and anything that comes before
        species_name=${species_name#*uncultured }
        # Removes first instance of ' bacterium' from the left side and anything that comes after
        species_name=${species_name%* bacterium}

        # Searching name.dmp first for a specific match, then for any matching species name and isolates the taxID of the first match
        taxid=$(awk -v species_name="$species_name" -F '|' '$2 == species_name {print $1; exit}' "$name_dump")
            if [ -z "$taxid" ]; then
                taxid=$(awk -v species_name="$species_name" -F '|' '$2 ~ species_name {print $1; exit}' "$name_dump")
                [ -z "$taxid" ] || echo -e "Inexact match found for:\t$original_species_name"
            else
                echo -e "Exact match found for:\t$original_species_name"
            fi

        # Checks whether taxID has been found and if it has creates an entry for the species in gtdb.accession2taxid
        if [ -n "$taxid" ]; then
            # Gets accession version and ID for current entry
            accession_version="$(awk -v linenum="$i" 'NR == linenum {print $1}' "$headerfile")"
            accession=${accession_version%.*}
            # Creates gi number for entry
            gi=$((600000000 + entrynum))
            ((entrynum++))
            
            echo -e "$accession\t$accession_version\t$taxid\t$gi" >> "$gtdb_accession2taxid_file"
        else
            echo -e "No match found for:\t$original_species_name"
        fi
    done
    duration=$SECONDS
    echo -e "Process took:\t$((duration / 3600)) hours $((duration / 60 % 60)) minutes $((duration % 60)) seconds"
    
    # Count number of lines in gtdb.accession2taxid file
    new_entries=$(wc -l <"$gtdb_accession2taxid_file")
    # Compress gtdb.accession2taxid file
    gzip "$gtdb_accession2taxid_file"
    
    echo -e "gtdb.accession2taxid.gz can be found here:\t'$gtdb_accession2taxid_file'.gz\nEntries added to gtdb.accession2taxid.gz:\t$((new_entries - 1))"
    end_date="$(date +"%Y-%m-%d")"
    end_time=$(date +"%T")
    echo -e "\nEnd date:\t$end_date\nEnd time:\t$end_time"
fi

exit 0