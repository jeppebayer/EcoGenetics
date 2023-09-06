#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

touch /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_filelist.txt
tempfilelist=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_tempfilelist.txt

for vcf in /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_*.vcf; do
    echo "$vcf" >> "$tempfilelist"
done

sort -n "$tempfilelist" > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_filelist.txt
rm -f "$tempfilelist"

bcftools concat \
-f /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/temp/Ocin_NYS-F_filelist.txt \
-o /faststorage/project/EcoGenetics/people/Jeppe_Bayer/data/Ocin_NYS-F.vcf \
-O v

exit 0