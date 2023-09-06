#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

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

stdout_dir="$sampledir"/out

# Removes empty .out files
for file in "$stdout_dir"/*; do
    if [ ! -s "$file" ]; then
        rm -f "$file"
    fi
done

logfile="$sampledir"/log_"$namebase"_"$(date +"%Y-%m-%d")".txt
touch "$logfile"

# Tries to identify .out files which report errors
grep -qxF "!!!ATTENTION!!!" "$logfile" || echo "!!!ATTENTION!!!" >> "$logfile"

for file in "$stdout_dir"/*; do
    if grep -qi "ERROR" "$file" || grep -qi "EXCEPTION" "$file" || grep -qi "WARNING" "$file"; then
        grep -qxF "Possible error in $namebase" || echo -e "\nPossible error in $namebase" >> "$logfile"
        errorfile=$(basename "$file")
        id=${errorfile##*_*_*-}; id=${id%%.out}
        echo -e "\t$errorfile\tJob ID: $id" >> "$logfile"
    fi
done

[ "$(tail -n 1 "$logfile")" == "!!!ATTENTION!!!" ] && sed -i '$d' "$logfile"

exit 0