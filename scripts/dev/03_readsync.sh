#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

awk '{print $4}' /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test1.txt > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/synccol4.txt

IFS=":"

while read -ra line; do
    
    if [ "${line[4]}" == "0" ] && [ "${line[5]}" == "0" ]; then
        
        n=0

        for (( i=0; i<=3; i++ )); do

            [ ! "${line[$i]}" == "0" ] && (( n++ ))
        
        done

        if [ $((n)) -le 2 ]; then

            echo -e "${line[0]}\t${line[1]}\t${line[2]}\t${line[3]}"

        fi

    fi

done < /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/synccol4.txt > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/sync.txt

rm -f /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/synccol4.txt
