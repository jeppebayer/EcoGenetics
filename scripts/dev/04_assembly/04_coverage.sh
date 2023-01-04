#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 06:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/coverage-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary

fi

target="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Entomobrya_nicoleti/depth/Enic_alignment.sam"
WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly"
fullname=$(basename "$target")
name=${fullname%.*}

samtools depth \
-@ 9 \
-a \
-H \
-o "$WD"/"$(basename "$(dirname "$(dirname "$target")")")"/depth/"$name".depth \
"$target"

exit 0