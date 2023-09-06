#!/bin/bash

reference_genome="$1"
repbaserun="$2"
repmodrun="$3"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate repeatmasking
fi

# Combine full RepeatMasker result files
cat \
"$repbaserun"/"$(basename "$reference_genome")".cat.gz \
"$repmodrun"/"$(basename "$reference_genome")".cat.gz \
> "$(basename "$reference_genome")".full_mask.cat.gz

# Combine RepeatMasker tabular files
cat \
"$repbaserun"/"$(basename "$reference_genome")".out \
<(tail -n +4 "$repmodrun"/"$(basename "$reference_genome")".out) \
> "$(basename "$reference_genome")".full_mask.out

ProcessRepeats \
-lib /faststorage/project/EcoGenetics/BACKUP/database/giri_repbase/RepBaseCustom15.02.fasta/Arthropoda.rep \
"$(basename "$reference_genome")".full_mask.cat.gz

# Make GFF 3 of repeat regions from RepeatMasker output file
repeatmasker_out="$(basename "$reference_genome")".full_mask.out
filename="$(basename "$repeatmasker_out")"
filename=${filename%.*.*.*}

awk \
    -F " " \
    'BEGIN {OFS="\t"; print "##gff-version 3"}
    NR > 3
        {if ($9 == "C")
            {strand = "-"}
        else
            {strand = "+"}
        if ($12 ~ /\(/)
            {start = $14}
        else {start = $12}
        print $5, "RepeatMasker", "repeat_region", $6, $7, ".", strand, ".", "ID="$15";Name="$10";Class="$11";Family="$11";Target="$10" "start" "$13}' \
    "$repeatmasker_out" \
    > "$filename".repeats.gff

exit 0