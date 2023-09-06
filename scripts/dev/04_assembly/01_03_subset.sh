#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --time=30
#SBATCH --mem=36G
#SBATCH --cpus-per-task=6
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

dbtable=$(readlink -f "$1")
asmfile=$(readlink -f "$2")
group="$3"

name=${dbtable%.*}

rg -j 6 '^# ' "$dbtable" > "$name"_"$group".txt
rg -j 6 -F "${group^}" "$dbtable" >> "$name"_"$group".txt

awk '{OFS="\t"} {print $1}' "$name"_"$group".txt > "$(dirname "$dbtable")"/temp_"$group"_idlist.txt

count=
while read -r line; do
    if [ -z $count ]; then
        rg -j 6 -F -A 1 "$line" "$asmfile" > "$(dirname "$dbtable")"/"$group"_"$(basename "$asmfile")"
        count="1"
    else
        rg -j 6 -F -A 1 "$line" "$asmfile" >> "$(dirname "$dbtable")"/"$group"_"$(basename "$asmfile")"
    fi
done < "$(dirname "$dbtable")"/temp_"$group"_idlist.txt

rm -f "$(dirname "$dbtable")"/temp_"$group"_idlist.txt

exit 0