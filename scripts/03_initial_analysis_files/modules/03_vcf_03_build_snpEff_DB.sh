#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage(){
cat << EOF

Usage: 03_vcf_00_init.sh -m -r <FILE>

Make custom database entry for reference genome in snpEff database.
Requires reference genome and annotation file in GTF or GFF format.

    -r  FILE        Reference genome file in FASTA format (.fa | .fna)

OPTIONS:
    -x              Overwrite existing identical entry.
    -h              Show this message.

EOF
}

# ----------------- Configuration ----------------------------------------

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

# Whether or not to overwrite exiting entry
overwrite=false

# Sets extended globbing
shopt -s extglob

# Don't mind
send_to_queue=true

# ------------------ Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 'r:xqh' OPTION; do
    case "$OPTION" in
        
        r)
            if [ -e "$OPTARG" ] && [[ "$OPTARG" == *.@(fna|fa) ]]; then
                reference_genome="$(readlink -f "$OPTARG")"
                reference_genome_dir="$(dirname "$reference_genome")"
            else
                echo -e "\nDesignated reference genome is not in accepted FASTA format (.fa | .fna)\n"
                exit 80
            fi
            ;;
        x)
            overwrite=true
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

# ------------------ Main -------------------------------------------------

gtf_file=("$reference_genome_dir"/*.gtf)
if [ ! -e "${gtf_file[*]}" ]; then
    gff_file=("$reference_genome_dir"/*.@(gff|gff3|gff2))
    if [ ! -e "${gff_file[*]}" ]; then
        echo -e "\nCan't locate GTF nor GFF annotation file in $reference_genome_dir/\n"
        exit 90
    fi
fi

if [ "$send_to_queue" == true ]; then
    if [ "$overwrite" == true ]; then
        jid=$(sbatch \
            --parsable \
            --account=EcoGenetics \
            --chdir="$reference_genome_dir" \
            --time=60 \
            --mem=100G \
            --cpus-per-task=2 \
            --output="$reference_genome_dir"/build_snpEff_DB.log \
            "$(readlink -f "$0")" -r "$reference_genome" -x -q)
    else
        jid=$(sbatch \
            --parsable \
            --account=EcoGenetics \
            --chdir="$reference_genome_dir" \
            --time=60 \
            --mem=100G \
            --cpus-per-task=2 \
            --output="$reference_genome_dir"/build_snpEff_DB.log \
            "$(readlink -f "$0")" -r "$reference_genome" -q)
    fi
    echo -e "\nLog file can be found here: $reference_genome_dir/build_snpEff_DB.log\n"
elif [ "$send_to_queue" == false ]; then
    echo -e "----- Job ID is $SLURM_JOB_ID -----\n"
    gtf_file=("$reference_genome_dir"/*.gtf)
    if [ ! -e "${gtf_file[*]}" ]; then
        gff_file=("$reference_genome_dir"/*.@(gff|gff3|gff2))
        if [ ! -e "${gff_file[*]}" ]; then
            echo -e "\nCan't locate GTF nor GFF annotation file in $reference_genome_dir/\n"
            exit 90
        fi
        gff_file_string="${gff_file[*]}"
        gff_file_string="${gff_file_string%.*}"

        agat_convert_sp_gff2gtf.pl --gff "${gff_file[*]}" --gtf_version 2.2 -o "$gff_file_string".gtf
        
        gtf_file=("$gff_file_string".gtf)
    fi

    reference_genome_name=$(basename "${reference_genome[*]}")
    reference_genome_version=${reference_genome_name%.*}
    reference_genome_version=${reference_genome_version%_genomic*}
    species_name=$(basename "$reference_genome_dir")
    species_name_space=${species_name/_/ }

    # snpEff 
    snpeff_path=("$(dirname "$(dirname "$(which python)")")"/share/snpeff*)

    snpeffdata="${snpeff_path[*]}"/data
    # Creates specified directory if it doesn't exist
    [ -d "$snpeffdata" ] || mkdir "$snpeffdata"

    # Creates backup of original config file if none exists
    snpeffconfig_backup="${snpeff_path[*]}"/snpEff.config.backup
    if [ ! -e "$snpeffconfig_backup" ]; then
        cp "${snpeff_path[*]}"/snpEff.config "$snpeffconfig_backup"
    fi

    # Checks whether the would be database entry already exists
    db_duplicate=$(rg -xFn -m 1 "# $species_name_space genome, version $reference_genome_version" "${snpeff_path[*]}"/snpEff.config)

    if [ -n "$db_duplicate" ]; then
        if [ "$overwrite" == true ]; then
            db_duplicate_line="${db_duplicate%:*}"
            sed -i ""$((db_duplicate_line - 1))","$((db_duplicate_line + 1))"d" "${snpeff_path[*]}"/snpEff.config
        else
            echo -e "\nDatabase for '$species_name_space genome, version $reference_genome_version' seems to already exist. Use '-f' to overwrite\n"
            exit 1
        fi
    fi

    # Adds new entry to config file
    echo -e "\n# $species_name_space genome, version $reference_genome_version\n$reference_genome_version.genome : $species_name_space" >> "${snpeff_path[*]}"/snpEff.config

    speciesdata="$snpeffdata"/"$reference_genome_version"
    # Creates specified directory if it doesn't exist
    [ -d "$speciesdata" ] || mkdir "$speciesdata"

    # Extract CDS sequences named using transcript ID 
    agat_sp_extract_sequences.pl -g "${gtf_file[*]}" -f "${reference_genome[*]}" -t cds -o "$speciesdata"/cds.fa

    # Extract protein sequences
    agat_sp_extract_sequences.pl -g "${gtf_file[*]}" -f "${reference_genome[*]}" -t cds -p -o "$speciesdata"/protein.fa

    # Creates copy of GTF file, genome reference sequence file, CDS file and protein file in speciesdata directory
    cp "${gtf_file[*]}" "$speciesdata"/genes.gtf
    cp "${reference_genome[*]}" "$speciesdata"/sequences.fa

    # Using explicit execution of snpEff.jar to change default java vm setting for heap size
    java -Xms80G -Xmx80G -jar "${snpeff_path[*]}"/snpEff.jar build -gtf22 -v "$reference_genome_version"
    # snpEff build -gtf22 -v "$reference_genome_version"

    # List of custom added genome to DB
    customdbs="${snpeff_path[*]}"/customDBs.txt
    if [ ! -f "$customdbs" ]; then
        echo -e "#Organism_name\tGenome_name" > "$customdbs"
    fi

    # Checks whether the wouldbe database entry already exists
    db_duplicate=$(rg -xFn -m 1 "$species_name\t$reference_genome_version" "$customdbs")

    if [ -n "$db_duplicate" ]; then
        db_duplicate_line="${db_duplicate%:*}"
        sed -i "${db_duplicate_line}d" "$customdbs"
    fi

    # Adds custom DB to list
    echo -e "$species_name\t$reference_genome_version" >> "$customdbs"

    if [ -e "$gff_file_string".agat.log ]; then
        rm "$gff_file_string".agat.log
    fi
fi

exit 0