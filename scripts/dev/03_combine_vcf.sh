#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

touch /home/jepe/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_filelist.txt
filelist=/home/jepe/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_filelist.txt

for vcf in /home/jepe/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_*.vcf; do
    echo "$(basename "$vcf")" >> "$filelist"
done

bcftools concat \
-f "$filelist" \
-o /home/jepe/EcoGenetics/people/Jeppe_Bayer/data/Ocin_NYS-F.vcf \
-O v

exit 0