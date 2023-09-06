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

filelist="$temp"/"$namebase"_spectrum_filelist.txt
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
        echo "$temp"/"$namebase"_"$num".spectrum >> "$filelist"

    else

        num="$i"
        echo "$temp"/"$namebase"_"$num".spectrum >> "$filelist"
        
    fi
done

sort -n -o "$filelist" "$filelist"

files=$(cat "$filelist")

# Creates array containing 0's with a length matching half the number of alleles
alleles=$((A / 2))

sfsarray=()

for (( i=1; i<=alleles; i++ )); do

    sfsarray[i]=0

done

oldIFS=$IFS

for file in $files; do

    IFS=$'\t' read -d '' -r -a arraypart <"$file"
    for (( i=1; i<=alleles; i++ )); do
        j=$(( i - 1 ))
        sfsarray[i]=$((sfsarray[i] + arraypart[j]))
    done

    if [ "$keep_temp" == "N" ]; then
        rm -f "$file"
    fi 
    
    IFS=$oldIFS
done

rm -f "$filelist"

# Creates spectrum file
spectrumfile="${S%.*}".spectrum
touch "$spectrumfile"

for (( i=1; i<=alleles; i++ ));do

    if [ $i == $alleles ]; then
        echo "${sfsarray[$i]}" >> "$spectrumfile"   
    else
        echo -ne "${sfsarray[$i]}\t" >> "$spectrumfile"
    fi

done

exit 0