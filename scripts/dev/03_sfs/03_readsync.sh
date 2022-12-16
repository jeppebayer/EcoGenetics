#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

# Replaces ':' in file with tab
sed 's/:/\t/g' /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs/Ocin_NYS-F.sync > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/temp.sync

# Reads observations from temporary file
while read -ra line; do
    
    # Checks that position doesn't contain any observations in column eight and nine, indicating unknown base or deletion
    if [ "${line[7]}" == "0" ] && [ "${line[8]}" == "0" ]; then
        
        # Counts how many positions contain observations
        n=0
        for (( i=3; i<=6; i++ )); do

            [ ! "${line[$i]}" == "0" ] && (( n++ ))
        
        done

        # Proceeds only with reads containing max 2 unique observations and where none has only a single observation
        if [ $((n)) -le 2 ] && [[ ! "${line[3]}" == "1" && ! "${line[4]}" == "1" && ! "${line[5]}" == "1" && ! "${line[6]}" == "1" ]]; then

            echo -e "${line[0]}\t${line[1]}\t${line[2]}\t${line[3]}\t${line[4]}\t${line[5]}\t${line[6]}"

        fi

    fi

done < /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/temp.sync > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs/Ocin_NYS-F_filtered.sync

rm -f /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/temp.sync

exit 0