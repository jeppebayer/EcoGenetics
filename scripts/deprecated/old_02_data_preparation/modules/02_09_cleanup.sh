#!/bin/bash

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
dataprep=$6 # Path to script location
algo=$7 # Chosen algorithm

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

stdout_dir="$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/stdout_sbatch

# Removes empty .out files
for file in "$stdout_dir"/*; do
    if [ ! -s "$file" ]; then
        rm -f "$file"
    fi
done

logfile="$WD"/"$dataprep"/"$(basename "$SD")"/log_"$(basename "$SD")"_"$(date +"%Y-%m-%d")".txt

# Tries to identify .out files which report errors
grep -qxF "!!!ATTENTION!!!" "$logfile" || echo "!!!ATTENTION!!!" >> "$logfile"

for file in "$stdout_dir"/*; do
    if grep -qi "ERROR" "$file" || grep -qi "EXCEPTION" "$file" || grep -qi "WARNING" "$file"; then
        grep -qxF "Possible error in $(basename "$sample")" || echo -e "\nPossible error in $(basename "$sample")" >> "$logfile"
        errorfile=$(basename "$file")
        id=${errorfile##*_*_*-}; id=${id%%.out}
        echo -e "\t$errorfile\tJob ID: $id" >> "$logfile"
    fi
done

[ "$(tail -n 1 "$logfile")" == "!!!ATTENTION!!!" ] && sed -i '$d' "$logfile"

exit 0