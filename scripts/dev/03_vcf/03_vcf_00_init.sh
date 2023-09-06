#!/bin/bash

# ------------------ Description ------------------------------------------

# Pipeline to create VCF files from large sets of pooled re-sequencing data
# utilizing freebayes

# ----------------- Configuration ----------------------------------------

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command={print $2}')
# If run on the frontend:
else
    script=$(realpath "$0")
fi
script_path="$(dirname "$script")"

# Default working directory
WD="$(readlink -f "$PWD")"

# Creates list for the execution of jobs
num_of_jobs=2
for (( i=1; i<="$num_of_jobs"; i++ )); do
    joblist["$i"]=true
done

jobgraph(){
cat << EOF

          1. Create VCF
                | |
    2. Concatenate all VCF parts

EOF
}

# ------------------ Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 'f:d:j:h' OPTION; do
    case "$OPTION" in
        f)
            if [ -f "$OPTARG" ] && [[ "$OPTARG" == *.bam ]]; then
                targetbam=$(readlink -f "$OPTARG")
            else
                echo -e "\n$OPTARG is not a BAM file\n"
                exit 1
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
if [ -z "$targetbam" ]; then
    echo -e "\nNo BAM file has been provided\n"
    exit 3
fi

# Removes "/" from the end of directory path
if [ "${WD:-1}" == "/" ]; then
    length=${#WD}
    WD=${WD::length - 1}
fi

# ----------------- Functions --------------------------------------------

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

job1(){
    jobidlist[1]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --array=1-"$ncontig" \
                --time=480 \
                --mem=50G \
                --cpus-per-task=1 \
                --dependency=afterany:"$(dependency "")" \
                --output="$out"/01_create_vcf-%a-%j.out \
                "$script_path"/modules/03_vcf_01_create_vcf.sh "$targetbam" "$RG" "$temp" "$sample_name")
}

job2(){
    jobidlist[2]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$WD" \
                --time=240 \
                --mem=36G \
                --cpus-per-task=6 \
                --dependency=afterany:"$(dependency "1")" \
                --output="$out"/02_concatenate-%j.out \
                "$script_path"/modules/03_vcf_02_concatenate.sh "$WD" "$temp" "$sample_name")
}

# ------------------ Main -------------------------------------------------

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

sample_name=$(basename "$(dirname "$targetbam")")
species_name=$(basename "$(dirname "$(dirname "$targetbam")")")

# Finds references genome for designated species or exits if it cannot be found
shopt -s extglob
for ref in /faststorage/project/EcoGenetics/BACKUP/reference_genomes/"$species_name"/*.@(fna|fa|fasta); do
    if [ ! -e "$ref" ]; then
		echo -e "\nCannot locate reference genome in /faststorage/project/EcoGenetics/BACKUP/reference_genomes/$species_name/, format: .fna, .fa, .fasta\n"
		exit 4
    else
		RG=$ref
		break
    fi
done

temp="$WD"/temp
# Creates specified directory if it doesn't exist
[ -d "$temp" ] || mkdir -m 775 "$temp"

WD="$WD"/03_initial_analysis_files
# Creates specified directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD"/"$species_name"
# Creates specified directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD"/"$sample_name"
# Creates specified directory if it doesn't exist
[ -d "$WD" ] || mkdir -m 775 "$WD"

out="$WD"/out
# Creates specified directory if it doesn't exist
[ -d "$out" ] || mkdir -m 775 "$out"

# Holds list of all contig names from bam file
SQlist=$(samtools view -H "$targetbam" | rg '@SQ' | awk '{print $2}')

# Turns SQlist into an array to easilier access individual contig names
x=1
for line in $SQlist; do
    SQarray["$x"]="$line"
    ((x++))
done

# Number of contigs in reference
ncontig=${#SQarray[@]}

# Executes job if indicated in joblist
[ "${joblist[1]}" = true ] && job1

# Executes job if indicated in joblist
[ "${joblist[2]}" = true ] && job2

exit 0