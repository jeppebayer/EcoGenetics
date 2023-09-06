#!/bin/bash

# ----------------- Description ------------------------------------------

# Script for navigating and initializing various utilities and pipelines
# used during genome assembly

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 04_assembly_init.sh [MODULE] [OPTIONS]

This script is used for navigating and initializing various utilities and
pipelines used during genome assembly

MODULE (must be assigned):
    hifiadapterfilt | filt      Removes any remaining adapters from HiFi reads.
    kmer_analysis | kmer        Runs kmer analysis on specified genome file in
                                fastx format utilizing jellyfish and GenomeScope.
    hifiasm                     Creates hifiasm genome assembly from specified
                                PacBio  HiFi reads.
    busco                       Runs BUSCO analysis on specified genome file in
                                fastx format.
    purge_dups | purge          Removes haplotigs and overlaps from genome assembly
                                based on read depth.
    blobtools                   Using coverage and sequence blasting (NCBI BLAST+
                                and DIAMOND) to separate possible contaminants
                                from genome assembly.
    gfa_to_fasta | fasta        Creates a fasta format version of a .gfa formatted file.

OPTIONS:
    -h | --help                 Show this message.
    -v                          Show current version and recent version history.
    -u  STRING                  Your username. Used to automatically change environment.
                                Be sure to have the necessary environmnets installed.

EOF
}

# ----------------- Configuration ----------------------------------------

# Current version and rescent version history
version()
{
cat << EOF

Current Version: 1.0.2
Updated 20/02-2023

Version history:
    Version             Changelog
    1.0.0               Added version information
    1.0.1               Added custom timer for BUSCO analysis.
                        Uses a logaritmic function based on filesize.
    1.0.2               Improved file name conflict solving in HiFiAdapterFilt.

EOF
}

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

script_path="$(dirname "$script_path")"

if [ ! "$USER" == "jepe" ]; then
    echo -e "\nThe scripts used in this collection have a large series of dependencies."
    echo "To avoid conflicts or environments of too great a size, these dependencies have been split on three different environments."
    echo "The files to install these environments can be found here: $script_path/envs/."
    echo "Each environment is installed using the following command:"
    echo "'conda env create -f <environment_name>.yml'"
    echo -e "Then pass the option '-u <your_username>' along with the chosen module and related options\n"
fi

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
    hifiadapterfilt | filt)
        bash "$script_path"/hifiadapterfilt/04_hifiadapterfilt_00_init.sh "${options[@]}"
        ;;
    kmer_analysis | kmer)
        bash "$script_path"/kmer_analysis/04_kmer_00_init.sh "${options[@]}"
        ;;
    hifiasm)
        bash "$script_path"/hifiasm/04_hifiasm_00_init.sh "${options[@]}"
        ;;
    busco)
        bash "$script_path"/busco/04_busco_00_init.sh "${options[@]}"
        ;;
    purge_dups | purge)
        bash "$script_path"/purge_dups/04_purge_dups_00_init.sh "${options[@]}"
        ;;
    blobtools)
        bash "$script_path"/blobtools/04_blobtools_00_init.sh "${options[@]}"
        ;;
    gfa_to_fasta | fasta)
        bash "$script_path"/gfa_to_fasta/04_gfa_to_fasta_00_init.sh "${options[@]}"
        ;;
    -h | --help)
        usage
        exit 0
        ;;
    -v)
        version
        exit 0
        ;;
    *)
        usage
        exit 1
        ;;
esac

exit 0