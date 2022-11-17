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

qualimap_path=$(dirname "$(which python)")/qualimap

line=$(awk '{if(NR==47) print $0}' "$qualimap_path")

echo "$line"

grep -qxF '\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE -XX:MaxPermSize=1024m"' "$qualimap_path" || sed -i '47s#.*#\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE -XX:MaxPermSize=1024m"#' "$qualimap_path"

line=$(awk '{if(NR==47) print $0}' "$qualimap_path")

echo "$line"

exit 0