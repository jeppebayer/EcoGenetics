#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

syncfile="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs/Ocin_NYS-F_filtered.sync"
A=100

# Creates array containing 0's with a length matching half the number of alleles
alleles=$((A / 2))

sfsarray=()

for (( i=1; i<=alleles; i++ )); do

    sfsarray[i]=0

done

# Reads each line from the filtered sync file one line at a time
while read -ra line; do

    # Assign the two non-zero observations to num1 and num2
    num1=

    for (( i=3; i<=6; i++ )); do

        if [ ! "${line[$i]}" == "0" ] && [ "$num1" ]; then

            num2=${line[$i]}
            break

        elif [ ! "${line[$i]}" == "0" ]; then

            num1=${line[$i]}

        fi

    done
        
    # Calculates total number of observations
    total=$(( num1 + num2 ))

    # Calculates proportion for each observation out of the total
    pct1=$(awk -v i="$num1" -v j="$total" 'BEGIN { print ( i / ( j / 100 ) ) }')
    pct2=$(awk -v i="$num2" -v j="$total" 'BEGIN { print ( i / ( j / 100 ) ) }')
    
    # Finds Minor Alleles count and round to nearest integer
    MA=$(awk -v i="$pct1" -v j="$pct2" -v k="$total" 'BEGIN { printf "%.0f", ( i < j )? i : j }')

    # Adds one to count for relevant Minor Allele count
    ((sfsarray[MA]+=1))

done < $syncfile

# Creates and writes to spectrum file
spectrum="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs/Ocin_NYS-F.spectrum"
touch "$spectrum"

for (( i=1; i<=alleles; i++ ));do

    if [ $i == $alleles ];then

        echo "${sfsarray[$i]}" >> $spectrum
    
    else

        echo -ne "${sfsarray[$i]}\t" >> $spectrum

    fi

done

# Removes empty .out files
for file in /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_sfs/out/*; do
    if [ ! -s "$file" ]; then
        rm -f "$file"
    fi
done

exit 0