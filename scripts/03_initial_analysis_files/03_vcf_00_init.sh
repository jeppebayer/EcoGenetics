#!/bin/bash

# ------------------ Description ------------------------------------------

# Pipeline to create VCF files from large sets of pooled re-sequencing data
# utilizing freebayes

# ----------------- Usage ------------------------------------------------

usage(){
cat << EOF

Usage: 03_vcf_00_init.sh -f <FILE> | -l <FILE> [OPTIONS]

Creates VCF file from BAM file using Freebayes. Further annotates VCFs using
snpEff and custom databases.

    -f  FILE        BAM file.
    -l  FILE        File containing list of BAM files to be analyzed

OPTIONS:
    -d  PATH        Path to working directory. If none is provided the
                    current working directory will be used.
    -p  INTEGER     Set ploidy of for the analysis. When samples are pooled,
                    it should equal the number of genome copies in the sample.
                    Default [100]
    -n  INTEGER     Sets the 'use-best-n-alleles' parameter. Sets the number
                    of haplotypes to consider for a given position. Seems
                    to only have a real effects where indels are considered.
                    Defualt [3].
    -c              Show list of custom databases.
    -m              Make custom SnpEff database.
    -i              Ignore custom database check.
    -r  FILE        Set reference genome manually.
    -j  STRING      If 'list' shows jobgraph. Can be supplied with a single
                    integer, 'x', or a range of integers 'x-y' to define 
                    which jobs in the pipeline should be run.
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

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script=$(scontrol show job "$SLURM_JOB_ID" | awk -F='' '/Command=/{print $2}' | awk '{print $1}')
# If run on the frontend:
else
    script=$(realpath "$0")
fi
script_path="$(dirname "$script")"

# Default working directory
WD="$(readlink -f "$PWD")"

# Default ploidy value
ploidy=100

# Default 'use-best-n-alleles' value (0 means all)
bestn=3

# Default for whether or not to skip check of custom database
skip_customdb_check=false

# Creates list for the execution of jobs
num_of_jobs=5
for (( i=1; i<="$num_of_jobs"; i++ )); do
    joblist["$i"]=true
done

jobgraph(){
cat << EOF

        1. Create VCF parts                          ________________
                | |                                 |                |
    2. Concatenate all VCF parts                    | Make SnpEff DB |
                | |                                 |________________|
 3. Annotation of VCF using snpEff
                | |
        4. Create Folded SFS

EOF
}

# ------------------ Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

if [ "$#" -ge 2 ]; then
    options=("$@")
    options=("${options[@]:1}")
fi

while getopts 'f:l:d:j:p:n:cir:mh' OPTION; do
    case "$OPTION" in
        f)
            if [ -e "$OPTARG" ] && [[ "$OPTARG" == *.bam ]]; then
                targetbam=$(readlink -f "$OPTARG")
            else
                echo -e "\n$OPTARG is not a BAM file\n"
                exit 1
            fi
            ;;
        l)
            if [ -e "$OPTARG" ] && [[ "$OPTARG" == .txt ]]; then
                targetbamlist=$(readlink -f "$OPTARG")
            else
                echo -e "\nList of BAM filest must be .txt\n"
                exit 20
            fi
            ;;
        d)
            WD=$(readlink -f "$OPTARG")
            ;;
        j)
            if [ "$OPTARG" == "list" ]; then
                jobgraph
                exit 0
            else
                for (( i=1; i<="$num_of_jobs"; i++ )); do
                    joblist["$i"]=false
                done
                if [ ${#OPTARG} -gt 2 ]; then
                    min=${OPTARG%-*}
                    max=${OPTARG#*-}
                    for (( i="$min"; i<="$max"; i++ )); do
                        joblist["$i"]=true
                    done
                else
                    joblist["$OPTARG"]=true
                fi
            fi
            ;;
        p)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                ploidy="$OPTARG"
            else
                echo -e "\nPloidy can only be assigned an integer value\n"
            fi
            ;;
        n)
            if [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                bestn="$OPTARG"
            else
                echo -e "\n'Use-best-n-alleles' can only be assigned an integer value\n"
            fi
            ;;
        c)
            snpeff_path=("$(dirname "$(dirname "$(which python)")")"/share/snpeff*)
            customdbs="${snpeff_path[*]}"/customDBs.txt
            if [ -e "$customdbs" ]; then
                cat <(echo) "$customdbs" <(echo)
                exit 0
            else
                echo -e "\nThere seems to be no custom Databases\n"
                exit 4
            fi
            ;;
        m)
            bash "$script_path"/modules/03_vcf_03_build_snpEff_DB.sh "${options[@]}"
            exit 0
            ;;
        i)
            skip_customdb_check=true
            ;;
        r)
            if [ -e "$OPTARG" ]; then
                RG=$(readlink -f "$OPTARG")
            else
                echo -e "\nDesignated reference genome $RG does not seem to exist\n"
                exit 10
            fi
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

# Checks that BAM file has been provided
if [ -z "$targetbam" ] && [ -z "$targetbamlist" ]; then
    echo -e "\nNo BAM file or list of BAM files has been provided\n"
    exit 3
fi

# Removes "/" from the end of directory path
if [ "${WD: -1}" == "/" ]; then
    length=${#WD}
    WD=${WD::length - 1}
fi

# ----------------- Functions --------------------------------------------

# Adjusts time based on the size of the largest region
timer () {
    first_region_entries=$(samtools view -c "$targetbam" "${SQarray[1]:3}")
    adjusted_time=$(awk -v first_region_entries="$first_region_entries" -v base_time="$1" 'BEGIN { print int( ( ( first_region_entries / 2934916 ) * 1.20 ) * base_time ) } ')
    echo "$adjusted_time"
}

# Function to flexibly control the dependencies of each job
dependency(){
    depend_on=
    for id in $1; do
         if [ -n "${jobidlist[$id]}" ]; then
            if [ -z "$depend_on" ]; then
                depend_on="${jobidlist[$id]}"
            else
                depend_on="$depend_on:${jobidlist[$id]}"
            fi
        fi
    done
    [ -z "$depend_on" ] && depend_on="-1"
    echo "$depend_on"
}
# "$(timer 600)"
# --array=1-"$nSQ"
# job1(){
#     jobidlist[1]=$(sbatch \
#                 --parsable \
#                 --account=EcoGenetics \
#                 --chdir="$WD" \
#                 --array=1-"$nSQ" \
#                 --time=4320 \
#                 --mem=50G \
#                 --cpus-per-task=1 \
#                 --dependency=afterany:"$(dependency "")" \
#                 --output="$out"/01_create_vcf_from_single-%A-%a-%j.out \
#                 "$script_path"/modules/03_vcf_01_create_vcf_from_single.sh "$targetbam" "$RG" "$temp" "$sample_name" "$bestn" "$ploidy")
# }

job1(){
    jobidlist[1]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --array=1-"$nSQ" \
                --time=03-00:00:00 \
                --mem=50G \
                --cpus-per-task=1 \
                --dependency=afterany:"$(dependency "")" \
                --output="$out"/01_create_vcf_from_single-%A-%a-%j.out \
                "$script_path"/modules/03_vcf_01_create_vcf_from_single.sh "$targetbam" "$RG" "$temp" "$sample_name" "$bestn" "$ploidy")
}

job2(){
    jobidlist[2]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --array=1-"$nSQ" \
                --time=960 \
                --mem=80G \
                --cpus-per-task=1 \
                --dependency=afterany:"$(dependency "")" \
                --output="$out"/01_create_vcf_from_list-%A-%a-%j.out \
                "$script_path"/modules/03_vcf_01_create_vcf_from_list.sh "$targetbamlist" "$RG" "$temp" "$sample_name" "$bestn" "$ploidy")
}

job3(){
    jobidlist[3]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time=240 \
                --mem=48G \
                --cpus-per-task=8 \
                --dependency=afterany:"$(dependency "1 2")" \
                --output="$out"/02_concatenate-%j.out \
                "$script_path"/modules/03_vcf_02_concatenate.sh "$WD" "$temp" "$sample_name" "$bestn")
}

job4(){
    jobidlist[4]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time=120 \
                --mem=50G \
                --cpus-per-task=2 \
                --dependency=afterany:"$(dependency "3")" \
                --output="$out"/03_snpEff-%j.out \
                "$script_path"/modules/03_vcf_03_snpEff.sh "$WD" "$sample_name" "$bestn" "$dbgenome")
}

job5(){
    jobidlist[5]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time=40 \
                --mem=50G \
                --cpus-per-task=2 \
                --dependency=afterany:"$(dependency "4")" \
                --output="$out"/04_pooledfoldedsfs-%j.out \
                "$script_path"/modules/03_vcf_04_pooledfoldedsfs.sh "$WD" "$sample_name" "$bestn" "$ploidy")
}

# ------------------ Main -------------------------------------------------
if [ -n "$targetbam" ]; then
    sample_name=$(basename "$(dirname "$targetbam")")
    species_name=$(basename "$(dirname "$(dirname "$targetbam")")")
else
    sample_name=$(basename "$(dirname "$(head -n 1 "$targetbamlist")")")
    # This shit is stupid... should overhaul entire naming convention
    species_name=$(basename "$(dirname "$(dirname "$(dirname "$(head -n 1 "$targetbamlist")")")")")
fi

# Finds references genome for designated species or exits if it cannot be found
shopt -s extglob
if [ -z "$RG" ]; then
    for ref in /faststorage/project/EcoGenetics/BACKUP/reference_genomes/"$species_name"/*.@(fna|fa|fasta); do
        if [ ! -e "$ref" ]; then
            echo -e "\nCannot locate reference genome in /faststorage/project/EcoGenetics/BACKUP/reference_genomes/$species_name/, format: .fna, .fa, .fasta\n"
            exit 5
        else
            RG=$ref
            break
        fi
    done
fi

if [ "$skip_customdb_check" = false ]; then
    # Checks whether a custom database file exists
    snpeff_path=("$(dirname "$(dirname "$(which python)")")"/share/snpeff*)
    customdbs="${snpeff_path[*]}"/customDBs.txt
    if [ ! -e "$customdbs" ]; then
        echo -e "\nThere semms to be no custom databases.\nWill abort script.\n"
        exit 4
    fi

    # Checks if an entry in the custom database file matches the reference organism name
    reference_genome_name="$(basename "$(dirname "$RG")")"
    dbgenome=$(awk -v reference_genome_name="$reference_genome_name" '$1 == reference_genome_name {print $2}' "$customdbs")
    if [ -z "$dbgenome" ]; then
        echo -e "\nThere seems to be no custom database entry for $reference_genome_name.\nWill abort script.\n"
        exit 6
    fi
fi

# temp="$WD"/temp
# # Creates specified directory if it doesn't exist ----------------------------------------------------- CHECK
# [ -d "$temp" ] || mkdir -m 775 "$temp"

WD="$WD"/03_initial_analysis_files
# Creates specified directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD"/"$species_name"
# Creates specified directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD"/"$sample_name"
# Creates specified directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

# Temporarily create temp directory within sample directory ---------------------------------------------- CHECK
temp="$WD"/temp
# Creates specified directory if it doesn't exist
[ -d "$temp" ] || mkdir -m 775 "$temp"

out="$WD"/out
# Creates specified directory if it doesn't exist
[ -d "$out" ] || mkdir -m 775 "$out"

# Holds list of all SQ names and LNs from bam file
if [ -n "$targetbam" ]; then
    SQlist=$(samtools view -H "$targetbam" | rg '@SQ' | awk '{print $2}')
else
    first_bam=$(head -n 1 "$targetbamlist")
    SQlist=$(samtools view -H "$first_bam" | rg '@SQ' | awk '{print $2}')
fi

# Turns SQlist and LNlist into an arrays to easilier access individual SQ names and LNs
x=1
for line in $SQlist; do
    SQarray["$x"]="$line"
    ((x++))
done

# Number of SQ in reference
nSQ=${#SQarray[@]}

# bpsperround=100000

if [ -n "$targetbam" ]; then
    # Executes job if indicated in joblist
    [ "${joblist[1]}" = true ] && job1
else
        # Executes job if indicated in joblist
    [ "${joblist[2]}" = true ] && job2
fi

# Executes job if indicated in joblist
[ "${joblist[3]}" = true ] && job3

# Executes job if indicated in joblist
[ "${joblist[4]}" = true ] && job4

# Executes job if indicated in joblist
[ "${joblist[5]}" = true ] && job5


exit 0