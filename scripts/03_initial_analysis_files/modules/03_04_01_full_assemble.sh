#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 30:00:00
#SBATCH --output=/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps/03_initial_analysis_files/03_sfs_full-%j.out

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
data=$6 # Path to data location in WD
n=$7 # Number of parts
script_path=$8 # Path to script location

# IFS=$'\t' read -r -a part < /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test1.txt

# echo "${part[0]}"

file1=/home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test1.txt
file2=/home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test2.txt

for part in {1..2..1}; do
    length=${#part}
    if [ "$length" -lt 2 ]; then
        num="$part"
        if [ "$num" == "1" ]; then
            cat /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test"$num".txt > /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test3.txt
        else
            cat /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test"$num".txt >> /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test3.txt
        fi
    else
        num="$part"

    fi
done

for txt in /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test*.txt; do
    
    while IFS=$'\t' read -r -a myarray; do

        for element in ${#myarray[@]}; do
            newarray+=("$element")
        done

    done < "$txt"

done
echo "${newarray[@]}"

exit 0