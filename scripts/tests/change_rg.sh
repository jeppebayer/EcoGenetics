#!/bin/bash

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate data_prep
fi

datapath=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/X_evolution/data/individual_data

if [ -z "$2" ]; then
    for bamfile in "$datapath"/*/*/*.bam; do
        sbatch \
            --parsable \
            --account=EcoGenetics \
            --chdir="$(dirname "$bamfile")" \
            --time=05:00:00 \
            --mem=80G \
            --cpus-per-task=10 \
            --output=/dev/null \
            --error=/dev/null \
            "$(readlink -f "$0")" "$bamfile" "a"
    done
else
    targetbam="$1"
    bamname="$(basename "$targetbam")"
    id="${bamname%_*.bam}"
    path="$(dirname "$targetbam")"
    bamnamenoext=${bamname%*.bam}

    samtools addreplacerg \
        -@ 9 \
        -m overwrite_all \
        -r "@RG\tID:$id\tSM:$id" \
        -O BAM \
        -o "$path"/"$bamnamenoext"_retag.bam \
        "$targetbam"

    samtools index \
        -b \
        -@ 9 \
        -o "$path"/"$bamnamenoext"_retag.bam.bai \
        "$path"/"$bamnamenoext"_retag.bam
fi

exit 0