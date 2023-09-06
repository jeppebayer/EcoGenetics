#!/bin/bash

sample="$1" # Sample directory
out="$2" # Out directory
logfile="$3" # Logfile
starttime="$4" # Pipeline starttime

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

# Removes empty .out files
for file in "$out"/*; do
    if [ ! -s "$file" ]; then
        rm -f "$file"
    fi
done

endtime="$(date +"%Y-%m-%d")"
if [ ! "$starttime" == "$endtime" ]; then
    echo -e "----- $endtime -----" >> "$logfile"
fi

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

[ "$(tail -n 1 "$logfile")" == "!!!ATTENTION!!!" ] && sed -i '$d' "$logfile" | sed -i '$d'

exit 0