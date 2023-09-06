#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --chdir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 2
#SBATCH --time 10:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/busco-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate busco

fi

asm_file="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Lepidocyrtus_sp/purge_dups/03/purged.fa"
speciesabbr="Lela"

busco \
-i "$asm_file" \
-m genome \
-o busco_"$speciesabbr" \
--out_path "$(dirname "$asm_file")" \
--download_path /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/busco/busco_downloads \
-l /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/busco/busco_downloads/lineages/arthropoda_odb10

exit 0