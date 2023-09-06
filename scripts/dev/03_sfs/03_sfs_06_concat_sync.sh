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

filelist="$temp"/"$namebase"_sync_filelist.txt
touch "$filelist"

# Adjusts naming according to file number
lineslength=${#lines}

for (( i=1; i<="$lines"; i++ )); do
    ilength=${#i}
    if [ $((ilength)) -lt $((lineslength)) ]; then

        dif=$((lineslength - ilength))
        num=""

        for (( j=1;  j<="$dif"; j++ )); do

            num="0$num"

        done

        num="$num$i"
        echo "$temp"/"$namebase"_"$num".sync >> "$filelist"

    else

        num="$i"
        echo "$temp"/"$namebase"_"$num".sync >> "$filelist"
        
    fi
done

sort -n -o "$filelist" "$filelist"

touch "$sampledir"/"$namebase".sync

files=$(cat "$filelist")

for file in $files; do

    cat "$file" >> "$sampledir"/"$namebase".sync

    if [ "$keep_temp" == "N" ]; then
        rm -f "$file"
        rm -f "$file".params
    fi 

done

rm -f "$filelist"

exit 0