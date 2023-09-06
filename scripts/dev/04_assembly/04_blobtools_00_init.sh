#!/bin/bash

asm_file="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Isotoma_sp/purge_dups/03/purged.fa"
hifi_reads="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Isotoma_sp/ISOSP_m64101e_221205_183109.hifi_reads.filt.fastq.gz"
species="Isotoma_sp"

WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps"
temp="$WD/temp"

WD="$WD/04_assembly"
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD/$species"
[ -d "$WD" ] || mkdir -m 775 "$WD"

WD="$WD/blobtools_nonpurged"
[ -d "$WD" ] || mkdir -m 775 "$WD"

out="$WD/out"
[ -d "$out" ] || mkdir -m 775 "$out"

firstletter=${species:0:1}
lastword=${species#*_}
threeletters=${lastword:0:3}
speciesabbr="$firstletter$threeletters"

if [ "$species" == "Isotoma_sp" ]; then
    speciesabbr="Iso_sp"
elif [ "$species" == "Lepidocyrtus_sp" ]; then
    speciesabbr="Lep_sp"
fi

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
    --time=180 \
    --mem-per-cpu=15G \
    --cpus-per-task=16 \
    --output="$out"/align-%j.out \
    "$script_path"/04_blobtools_01_align.sh "$asm_file" "$hifi_reads" "$WD" "$temp" "$speciesabbr" )

jid2=$(sbatch \
    --parsable \
    --time=180 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid1" \
    --output="$out"/cov-%j.out \
    "$script_path"/04_blobtools_02_cov.sh "$asm_file" "$hifi_reads" "$WD" "$temp" "$speciesabbr" )

jid3_1=$(sbatch \
    --parsable \
    --time=4320 \
    --mem-per-cpu=15G \
    --cpus-per-task=16 \
    --output="$out"/blastn-%j.out \
    "$script_path"/04_blobtools_03_blastn.sh "$asm_file" "$hifi_reads" "$WD" "$temp" "$speciesabbr" )

jid3_2=$(sbatch \
    --parsable \
    --time=4320 \
    --mem-per-cpu=15G \
    --cpus-per-task=16 \
    --output="$out"/blastd-%j.out \
    "$script_path"/04_blobtools_03_blastd.sh "$asm_file" "$hifi_reads" "$WD" "$temp" "$speciesabbr" )

jid4=$(sbatch \
    --parsable \
    --time=600 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid2":"$jid3_1":"$jid3_2" \
    --output="$out"/create-%j.out \
    "$script_path"/04_blobtools_04_create.sh "$asm_file" "$hifi_reads" "$WD" "$temp" "$speciesabbr" )

# jid4=$(sbatch \
#     --parsable \
#     --time=600 \
#     --mem-per-cpu=10G \
#     --cpus-per-task=2 \
#     --dependency=afterany:"$jid3_1":"$jid3_2" \
#     --output="$out"/create-%j.out \
#     "$script_path"/04_blobtools_04_create.sh "$asm_file" "$hifi_reads" "$WD" "$temp" "$speciesabbr" )

jid5=$(sbatch \
    --parsable \
    --time=180 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid4" \
    --output="$out"/plot-%j.out \
    "$script_path"/04_blobtools_05_plot.sh "$asm_file" "$hifi_reads" "$WD" "$temp" "$speciesabbr" )

jid6=$(sbatch \
    --parsable \
    --time=120 \
    --mem-per-cpu=10G \
    --cpus-per-task=2 \
    --dependency=afterany:"$jid4" \
    --output="$out"/view-%j.out \
    "$script_path"/04_blobtools_06_view.sh "$asm_file" "$hifi_reads" "$WD" "$temp" "$speciesabbr" )

exit 0