#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

adjustment=0.1

# number=$(awk -v adjustment=$adjustment 'BEGIN { printf "%.0f\n", ( 720 * adjustment) }')
# number=$(awk -v adjustment=$adjustment 'BEGIN { print int( 720 * adjustment + 120) }')
# $((2/filesize))
# echo "$number"

# export DISPLAY=:0

# echo "$PWD"
# echo "$WD"

# qualimap_path=$(dirname "$(which python)")/qualimap

# line=$(awk '{if(NR==47) print $0}' "$qualimap_path")

# echo "$line"

# grep -qxF '\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE -XX:MaxPermSize=1024m"' "$qualimap_path" || sed -i '47s#.*#\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE -XX:MaxPermSize=1024m"#' "$qualimap_path"

# line=$(awk '{if(NR==47) print $0}' "$qualimap_path")

# echo "$line"

# mapfile -t id < <(squeue -u "$USER" -S -V -o %A-%V)

# echo "${id[1]}"

# parts=20

# parts=$(awk -v p="$parts" 'BEGIN { print ( p - 1 ) }')

# echo "$parts"

# last=$(tail -n 2 /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/tests/test2.sh); last="${last%%$'\n'*}"

# errorfile=AgUr_01_J_2022_07_06-10434619

# id=${errorfile##*_*_*-}; id=${id%%.out}

# echo "$id"

# errorfile=Ocin_NYS-F_2022_07_06-10434619

# id=${errorfile#*_*_*-}
# id=${id%%.out}

# echo "$id"

# WD=/home/jepe/EcoGenetics/people/Jeppe_Bayer/steps

# sample=/home/jepe/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/NYS-F_C25

# R1_1=
# for R1_2 in "$sample"/*_1.fq.gz ; do
#     if [ "$R1_1" ]; then
#         cat "$R1_1" "$R1_2" > "$WD"/temp/"$(basename "$sample")"_R1.fq.gz
#     else
#         R1_1=$R1_2
#     fi
# done

env_path="$(dirname "$(dirname "$(dirname "$(which python)")")")"

echo "$env_path/poolhmm"

# a="1"
# b=${#a}

# if [ $(($b)) -lt 2 ]; then
#     echo "0$a"
# else
#     echo "$a"
# fi

exit 0