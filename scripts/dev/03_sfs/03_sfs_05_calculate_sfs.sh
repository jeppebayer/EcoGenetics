#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

lines=$1 # Number of unique IDs
qname=$2 # File containing IDs
S=$3 # Path to sample .mpileup file
A=$4 # Number of alleles in sample
namebase=$5 # Current name base for sample
temp=$6 # Directory for temporary files
sampledir=$7 # Sample specific directory within data directory
keep_temp=$8 # Flag for whether or not temporary files should be kept

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary

fi

# Creates array containing 0's with a length matching half the number of alleles
alleles=$((A / 2))

sfsarray=()

for (( i=1; i<=alleles; i++ )); do

    sfsarray[i]=0

done

# Function to count minor alleles
count_minor()
{
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

    done < "$temp"/"$namebase"_filtered_"$num".sync
}

# Function to create and write to spectrum part
create_spectrum()
{
    spectrum="$temp"/"$namebase"_"$num".spectrum
    touch "$spectrum"

    for (( i=1; i<=alleles; i++ ));do

        if [ $i == $alleles ]; then

            echo "${sfsarray[$i]}" >> "$spectrum"
        
        else

            echo -ne "${sfsarray[$i]}\t" >> "$spectrum"

        fi

    done
}

# Adjusts naming according to file number
id=${SLURM_ARRAY_TASK_ID}
idlength=${#id}
lineslength=${#lines}

if [ $((idlength)) -lt $((lineslength)) ]; then

    dif=$((lineslength - idlength))
    num=""

    for (( i=1;  i<="$dif"; i++ )); do

        num="0$num"

    done

    num="$num$id"
    count_minor
    create_spectrum

else

    num="$id"
    count_minor
    create_spectrum
    
fi

exit 0