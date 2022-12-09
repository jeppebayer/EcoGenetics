#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 1
#SBATCH --time 01:00:00

lines=$(wc -l < /home/jepe/EcoGenetics/people/Jeppe_Bayer/data/qname.txt)

id="1"
idlength=${#id}
lineslength=${#lines}

if [ $((idlength)) -lt $((lineslength)) ]; then
    dif=$((lineslength - idlength))
    num=""
    for (( i=1;  i<="$dif"; i++ )); do
        num="0$num"
    done
    num="$num$id"
else
    num="$id"
fi

exit 0