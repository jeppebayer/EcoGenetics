#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 00:40:00
#SBATCH --output=SFS_Full-%j.out

temp="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp"
scripts="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/03_sfs"
WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs"
out="$WD"/out

[ -d "$temp" ] || mkdir -m 775 "$temp"
[ -d "$WD" ] || mkdir -m 775 "$WD"
[ -d "$out" ] || mkdir -m 775 "$out"

# Creates file with all QNAME from .bam file
samtools view -H /faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F/Ocin_NYS-F_filtered.bam \
| grep '@SQ' \
| awk '{print $2}' \
> "$temp"/headertemp.txt

while read -r line; do
    echo "${line:3}"
done < "$temp"/headertemp.txt > "$temp"/qname.txt

# Removes temp header file
rm -f "$temp"/headertemp.txt

# Number of lines, eg. number of unique IDs
lines=$(wc -l < "$temp"/qname.txt)

jid1=$(sbatch \
    --parsable \
    --array=1-"$lines" \
    --time=120 \
    --mem-per-cpu=10G \
    --cpus-per-task=5 \
    --output="$out"/split_pileup-%a-%j.out \
    "$scripts"/03_split_pileup.sh "$lines" "$temp"/qname.txt)

jid2=$(sbatch \
    --parsable \
    --array=1-"$lines" \
    --time=180 \
    --mem-per-cpu=10G \
    --cpus-per-task=5 \
    --dependency=aftercorr:"$jid1" \
    --output="$out"/filter_coverage-%a-%j.out \
    "$scripts"/03_filter_coverage.sh "$lines" "$temp"/qname.txt)

jid3=$(sbatch \
    --parsable \
    --array=1-"$lines" \
    --time=180 \
    --mem-per-cpu=10G \
    --cpus-per-task=5 \
    --dependency=aftercorr:"$jid2" \
    --output="$out"/mpileup_to_sync-%a-%j.out \
    "$scripts"/03_mpileup_to_sync.sh "$lines" "$temp"/qname.txt)

jid4=$(sbatch \
    --parsable \
    --time=180 \
    --mem-per-cpu=10G \
    --cpus-per-task=10 \
    --dependency=afterany:"$jid3" \
    --output="$out"/concat_sync-%j.out \
    "$scripts"/03_concat_sync.sh "$lines" "$temp"/qname.txt)

jid5=$(sbatch \
    --parsable \
    --time=180 \
    --mem-per-cpu=10G \
    --cpus-per-task=10 \
    --dependency=afterany:"$jid4" \
    --output="$out"/readsync-%j.out \
    "$scripts"/03_readsync.sh "$lines" "$temp"/qname.txt)

exit 0