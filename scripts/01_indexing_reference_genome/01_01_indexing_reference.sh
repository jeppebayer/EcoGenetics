#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH --cpus-per-task 8
#SBATCH --time 2:00:00
#SBATCH --output=Index_Reference_Genome-%j.out

# Species specific reference genome (FASTA format)
RG=

# Option to force run even if target sample directory contains .bam
force_overwrite="N"

# ----------------- Script Flag Processing -------------------------------

while getopts 'r:fh' OPTION; do
    case "$OPTION" in
        r)
            RG="$(readlink -f "$OPTARG")"
            ;;
        f)
            force_overwrite="Y"
            ;;
        h)
            usage
            exit 1
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done

# 
if [ ! "$RG" == "*.fna" ]; then
    echo "$RG is not a .fna file"
    exit 1
fi

# Exit with error message if reference genome is indexed
for index in "$(dirname "$RG")"/*.ann; do
    if [ -e "$index" ] && [ "$force_overwrite" == "N" ]; then
        echo "Designated reference genome is already indexed, $RG"
        exit 1
    fi
done

# ----------------- Indexing ---------------------------------------------

# Indexing reference genome for aligning
bwa index \
-p "${RG%.*}" \
"$RG" ; \

# Adding fai index to reference
samtools faidx \
-o "${RG%.*}".fai \
"$RG"

exit 0