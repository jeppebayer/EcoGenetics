#!/bin/bash

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
dataprep=$6 # Path to script location
algo=$7 # Chosen algorithm

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

samtools depth \
-@ 7 \
-a \
-o "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.depth \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam

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

# Creates histgram of read depth
python "$script_path"/02_08_07_plotdepth.py \
"$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.depth \
"$sample"

exit 0