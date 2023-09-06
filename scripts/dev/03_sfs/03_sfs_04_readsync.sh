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

# Function to replace ':' in .sync file with tab
change_to_tab()
{
    sed 's/:/\t/g' "$temp"/"$namebase"_"$num".sync > "$temp"/"$namebase"_temp_"$num".sync
}

# Function to seads observations from temporary file
readsync()
{
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

    done < "$temp"/"$namebase"_temp_"$num".sync > "$temp"/"$namebase"_filtered_"$num".sync
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
    change_to_tab
    readsync

else

    num="$id"
    change_to_tab
    readsync
    
fi

rm -f "$temp"/"$namebase"_temp_"$num".sync

exit 0