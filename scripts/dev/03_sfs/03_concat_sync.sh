#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

filelist="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_sync_filelist.txt"
touch "$filelist"

for sync in /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_*.sync; do
    
    echo "$sync" >> "$filelist"

done

sort -n -o "$filelist" "$filelist"

touch /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs/Ocin_NYS-F.sync

files=$(cat "$filelist")

for file in $files; do

    cat "$file" >> /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs/Ocin_NYS-F.sync
    rm -f "$file"

done

rm -f "$filelist"

exit 0