#!/bin/bash

targetbam="$1" # Target bam file

targetbam=/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/OrcCin_NYS-F/OrcCin_NYS-F_filtered.bam

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

# awk -F "\t" 'END {print NR}' \
# <(awk -F "\t" '$3 >= 200 {print $0}' \
# <(samtools depth \
# -@ 7 \
# -a \
# "$targetbam"))

samtools depth \
-@ 7 \
-a \
"$targetbam" \
| awk -F "\t" '
$3 >= 200 {
    print $0
}' \
| awk -F "\t" '
END {
    print NR
}'

exit 0