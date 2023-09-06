#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 00_toolbox.sh [PIPELINE | MODULE] [OPTIONS]

Welcome to the TOOLBOX.
This script is used for navigating and initializing various utilities and
pipelines.

PIPELINE:
    01_indexing_reference_genome | reference
    02_data_preparation | dataprep
    03_initial_analysis_files | initanalysis
    04_assembly | assembly

MODULE:
    genome_length
    plotgc
    subset
    name_asm
    bam_to_fq
    gtdb_accession2taxid

OPTIONS:
    -h | --help                 Show this message.

EOF
}

# ----------------- Configuration ----------------------------------------

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

script_path="$(dirname "$script_path")"
pipe_path="$(dirname "$script_path")"

# ----------------- Module Processing ------------------------------------

# Makes sure all associated modules are executable
for script in "$script_path"/*; do
    [ ! -x "$script" ] && chmod -R u+rwx "$script_path" && break
done

if [ "$#" -ge 2 ]; then
    options=("$@")
    options=("${options[@]:1}")
fi

if [ -z "$1" ]; then
    usage
    exit 1
fi

case "$1" in
    01_indexing_reference_genome | reference)
        bash "$pipe_path"/01_indexing_reference_genome/01_initialize_indexing_reference.sh "${options[@]}"
        ;;
    02_data_preparation | dataprep)
        bash "$pipe_path"/02_data_preparation/02_initialize_data_prep.sh "${options[@]}"
        ;;
    03_initial_analysis_files | initanalysis)
        bash "$pipe_path"/03_initial_analysis_files/03_initialize_init_analysis_files.sh "${options[@]}"
        ;;
    04_assembly | assembly)
        bash "$pipe_path"/04_assembly/04_assembly_init.sh "${options[@]}"
        ;;
    genome_length)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: genome_length <PATH>\n\nCreates file containing total length in bp and a file with sorted file\nwith all sequence names and their length of a genome file in FASTA format\n\nFILE\tPath to genome file in FASTA format\n"
        else
            bash "$script_path"/modules/00_genome_length.sh "${options[@]}"
        fi
        ;;
    plotgc)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: plotgc <FILE> <VAL1> <VAL2>\n\nPlots histogram of GC content from BlobtoolsDB\n\n\tFILE\t\tPath to BlobtoolsDB file in .txt format\n\n\tVAL1\t\tOptional x intercept value for vertical line\n\tVAL2\t\tOptional x intercept value for vertical line\n"
        else
            if [ "$USER" == "jepe" ]; then
            # shellcheck disable=1090
            source /home/"$USER"/.bashrc
            # shellcheck disable=1091
            source activate R
            fi
            Rscript "$script_path"/modules/00_blobtools_plotgc.r "${options[@]}"
        fi
        ;;
    subset)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: subset -t <FILE> -a <FILE> [OPTIONS]\n\nSubsets genome assembly file based on Blobtools database.\n\n\t-t\tFILE\t\t\tBlobtoolsDB file in table text format.\n\t-a\tFILE\t\t\tGenome assembly file to subset.\n\nOPTIONS (Max. 1 option at a time):\n\t-g\tSTRING\t\t\tGroup name to filter on ex. 'arthropoda'.\n\t-e\tSTRING\t\t\tExclude specific group.\n\t-m\tRANGE\t\t\tRange of GC content to filter on in as 'x-y'.\n\t-c\tRANGE | MIN\t\tCoverage to filter on either a range in format\n\t\t\t\t\t'x-y' or a min. value in format 'x'.\n\t-l\tRANGE | MIN\t\tSequence lenght to filter on either a range in\n\t\t\t\t\tformat 'x-y' or a min. value in format 'x'.\n"
        else
        # shellcheck disable=2034
        jid1=$(sbatch \
            --parsable \
            --account=EcoGenetics \
            --time=30 \
            --mem=24G \
            --cpus-per-task=6 \
            --output=/dev/null \
            --error=/dev/null \
            "$script_path"/modules/00_blobtools_subset.sh "${options[@]}")
        fi
        ;;
    name_asm)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: name_asm -f <FILE> [OPTIONS]\n\nCreates standardized named version of assembly in reference genome directory\n\n\t-f\tFILE\t\t\tAssembly file.\n\t-m\t\t\t\tIf assembly has not been repeat masked.\n\t-a\t\t\t\tIf assembly has not been annotated.\n"
        else
            bash "$script_path"/modules/00_name_asm.sh "${options[@]}"
        fi
        ;;
    bam_to_fq)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: bam_to_fq -s <DIR> -d <DIR> [OPTIONS]\n\nCreates FASTQ formatted version of BAM file.\n\n\t-s\tDIRECTORY\t\t\tSpecies specific directory containing sample folders.\n\t-d\tDIRECTORY\t\t\tTarget directory in which to place FASTQ files.\n"
        else
            bash "$script_path"/modules/00_bam_to_fastq.sh "${options[@]}"
        fi
        ;;
    gtdb_accession2taxid)
        bash "$script_path"/modules/00_gtdb_accession2taxid.sh "${options[@]}"
        ;;
    -h | --help)
        usage
        exit 0
        ;;
    *)
        usage
        exit 1
        ;;
esac

exit 0