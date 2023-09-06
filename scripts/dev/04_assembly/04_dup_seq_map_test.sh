#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 6
#SBATCH --time 01:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Entomobrya_nicoleti/purge_dups/02/dup_seq/out/dup_seq_map_test-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools

fi

target="$(readlink -f "$1")"

temp="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp"

minimap2 \
-x map-hifi \
-a \
-t 6 \
"$target" \
"$target" \
| samtools sort \
-@ 15 \
-o "$(dirname "$target")"/"$(basename "$target")".bam \
-O bam \
-T "$temp" \
-

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary

fi

samtools index \
-b \
-@ 5 \
"$(dirname "$target")"/"$(basename "$target")".bam

samtools depth \
-a \
-@ 5 \
-o "$(dirname "$target")"/"$(basename "$target")".depth \
"$(dirname "$target")"/"$(basename "$target")".bam

exit 0