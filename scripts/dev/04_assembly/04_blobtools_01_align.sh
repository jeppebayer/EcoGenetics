#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools

fi

asm_file=$1
hifi_reads=$2
WD=$3
temp=$4
speciesabbr=$5

minimap2 \
-x map-hifi \
-a \
-t 16 \
"$asm_file" \
"$hifi_reads" \
| samtools sort \
-@ 15 \
-o "$WD"/"$speciesabbr".bam \
-O bam \
-T "$temp" \
-

samtools index \
-b \
-@ 15 \
"$WD"/"$speciesabbr".bam 

exit 0