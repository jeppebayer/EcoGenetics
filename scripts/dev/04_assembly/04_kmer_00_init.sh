#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 2
#SBATCH --time 00:40:00
#SBATCH --output=Kmer_init-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

target="$(readlink -f "$1")"

temp="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp"
WD="$(dirname "$target")/kmer"
out="$WD"/out

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

# Holds path to script directory
script_path=$(dirname "$script_path")

[ -d "$temp" ] || mkdir -m 775 "$temp"
[ -d "$WD" ] || mkdir -m 775 "$WD"
[ -d "$out" ] || mkdir -m 775 "$out"

jid1=$(sbatch \
    --parsable \
    --array="$2" \
    --time=40 \
    --mem-per-cpu=15G \
    --cpus-per-task=32 \
    --output="$out"/kmer_genome_est-%a-%j.out \
    "$script_path"/04_kmer_01_count.sh "$target" "$WD")

jid2=$(sbatch \
    --parsable \
    --array="$2" \
    --time=30 \
    --mem-per-cpu=15G \
    --cpus-per-task=32 \
    --dependency=aftercorr:"$jid1" \
    --output="$out"/kmer_histogram-%a-%j.out \
    "$script_path"/04_kmer_02_histogram.sh "$target" "$WD")

jid3=$(sbatch \
    --parsable \
    --array="$2"\
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=aftercorr:"$jid2" \
    --output="$out"/genomescope-%a-%j.out \
    "$script_path"/04_kmer_03_genomescope.sh "$target" "$WD" "$3")

exit 0