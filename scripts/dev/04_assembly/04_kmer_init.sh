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

target="/faststorage/project/EcoGenetics/BACKUP/PacBio_HiFi/Entomobrya_nicoleti/ENIC21_m64101e_221211_025112.hifi_reads.fastq.gz"

temp="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp"
scripts="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/dev/04_assembly"
WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/$(basename "$(dirname "$target")")"
out="$WD"/out

[ -d "$temp" ] || mkdir -m 775 "$temp"
[ -d "$WD" ] || mkdir -m 775 "$WD"
[ -d "$out" ] || mkdir -m 775 "$out"

jid1=$(sbatch \
    --parsable \
    --array=5-39:2 \
    --time=40 \
    --mem-per-cpu=15G \
    --cpus-per-task=32 \
    --output="$out"/kmer_genome_est-%a-%j.out \
    "$scripts"/04_kmer_genome_size_estimation.sh "$target" "$WD")

jid2=$(sbatch \
    --parsable \
    --array=5-39:2 \
    --time=30 \
    --mem-per-cpu=15G \
    --cpus-per-task=32 \
    --dependency=aftercorr:"$jid1" \
    --output="$out"/kmer_histogram-%a-%j.out \
    "$scripts"/04_kmer_histogram.sh "$target" "$WD")

jid3=$(sbatch \
    --parsable \
    --array=5-39:2\
    --time=30 \
    --mem-per-cpu=10G \
    --cpus-per-task=4 \
    --dependency=aftercorr:"$jid2" \
    --output="$out"/kmer_plot_hold-%a-%j.out \
    "$scripts"/04_kmer_plot_hold.sh "$target" "$WD")

exit 0