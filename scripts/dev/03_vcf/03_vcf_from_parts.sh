#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 00:40:00
#SBATCH --output=VCF_from_parts-%j.out

# Creates file with all QNAME from .bam file
samtools view -H /faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F/Ocin_NYS-F_filtered.bam \
| grep '@SQ' \
| awk '{print $2}' \
> /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/headertemp.txt

while read -r line; do
    echo "${line:3}"
done < /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/headertemp.txt > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/qname.txt

# Removes temp header file
rm -f /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/headertemp.txt

# Number of lines, eg. number of unique IDs
lines=$(wc -l < /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/qname.txt)

# jid1=$(sbatch \
#     --parsable \
#     --array=1-"$lines" \
#     --time=180 \
#     --mem-per-cpu=10G \
#     --cpus-per-task=10 \
#     --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/out/split_bam-%a-%j.out \
#     /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/03_split_bam.sh "$lines" people/Jeppe_Bayer/data/qname.txt)

# jid2=$(sbatch \
#     --parsable \
#     --array=1-"$lines" \
#     --time=300 \
#     --mem-per-cpu=30G \
#     --cpus-per-task=5 \
#     --dependency=afterany:"$jid1" \
#     --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/out/freebayes-%a-%j.out \
#     /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/03_freebayes.sh "$lines" people/Jeppe_Bayer/data/qname.txt)

jid2=$(sbatch \
    --parsable \
    --array=1-"$lines" \
    --time=600 \
    --mem-per-cpu=30G \
    --cpus-per-task=10 \
    --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/out/freebayes-%a-%j.out \
    /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/03_freebayes.sh "$lines" people/Jeppe_Bayer/data/qname.txt)

# sbatch \
#     --parsable \
#     --time=360 \
#     --mem-per-cpu=10G \
#     --cpus-per-task=10 \
#     --dependency=afterany:"$jid2" \
#     --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/out/concatenate-%a-%j.out \
#     /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/03_combine_vcf.sh