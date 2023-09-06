#!/bin/bash

target="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Lepidocyrtus_sp/LELA_m64101e_221212_121647.hifi_reads.filt.fastq.gz"
ref="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Lepidocyrtus_sp/purge_dups/02/purged.fa"
round="03"

WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly"
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD/$(basename "$(dirname "$target")")"
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD/purge_dups"
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD/$round"
[ -d "$WD" ] || mkdir -m 775 "$WD"

out="$WD/out"
[ -d "$out" ] || mkdir -m 775 "$out"

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

jid1=$(sbatch \
    --parsable \
    --time=60 \
    --mem-per-cpu=10G \
    --cpus-per-task=15 \
    --output="$out"/minimap2_1-%j.out \
    "$script_path"/04_purge_dups_01_minimap2.sh "$target" "$ref" "$WD")

jid2=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid1" \
    --output="$out"/pbcstat-%j.out \
    "$script_path"/04_purge_dups_02_pbcstat.sh "$target" "$ref" "$WD")

jid3=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid2" \
    --output="$out"/calcuts-%j.out \
    "$script_path"/04_purge_dups_03_calcuts.sh "$target" "$ref" "$WD")

jid4=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid3" \
    --output="$out"/split-%j.out \
    "$script_path"/04_purge_dups_04_split.sh "$target" "$ref" "$WD")

jid5=$(sbatch \
    --parsable \
    --time=60 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid4" \
    --output="$out"/minimap2_2-%j.out \
    "$script_path"/04_purge_dups_05_minimap2.sh "$target" "$ref" "$WD")

jid6=$(sbatch \
    --parsable \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=15 \
    --dependency=afterany:"$jid5" \
    --output="$out"/purge-%j.out \
    "$script_path"/04_purge_dups_06_purge.sh "$target" "$ref" "$WD")

jid7=$(sbatch \
    --parsable \
    --chdir="$WD" \
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid6" \
    --output="$out"/get_seqs-%j.out \
    "$script_path"/04_purge_dups_07_get_seqs.sh "$target" "$ref" "$WD")

exit 0