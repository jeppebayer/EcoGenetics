#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 15
#SBATCH --time 06:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/minimap2-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

target="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Entomobrya_nicoleti/ENIC21_m64101e_221211_025112.hifi_reads.filt.fastq.gz"
ref="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Entomobrya_nicoleti/Enic.asm.bp.p_ctg.fasta"
WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly"
name=$(basename "$target")
firstletter=${name:0:1}
nextletters="${name:1:3}"
lowerletters="${nextletters,,}"

# SAM output
# minimap2 \
# -x map-hifi \
# -t 15 \
# -a \
# -o "$WD"/"$(basename "$(dirname "$target")")"/"$firstletter""$lowerletters"_alignment.sam \
# "$ref" \
# "$target"

# PAF output
minimap2 \
-x map-hifi \
-t 15 \
-o "$WD"/"$(basename "$(dirname "$target")")"/"$firstletter""$lowerletters"_alignment.paf \
"$ref" \
"$target"

exit 0