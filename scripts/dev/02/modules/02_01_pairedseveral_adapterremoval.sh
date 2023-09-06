#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

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

# Concatenation of R1 files
R1_1=
for R1_2 in "$sample"/*_1.fq.gz ; do
    if [ "$R1_1" ]; then
        cat "$R1_1" "$R1_2" > "$WD"/temp/"$(basename "$sample")"_R1.fq.gz
    else
        R1_1=$R1_2
    fi
done

# Concatenation of R2 files
R2_1=
for R2_2 in "$sample"/*_2.fq.gz ; do
    if [ "$R2_1" ]; then
        cat "$R2_1" "$R2_2" > "$WD"/temp/"$(basename "$sample")"_R2.fq.gz
    else
        R2_1=$R2_2
    fi
done

# AdapterRemoval using concatenated files
AdapterRemoval \
--threads 8 \
--file1 "$WD"/temp/"$(basename "$sample")"_R1.fq.gz \
--file2 "$WD"/temp/"$(basename "$sample")"_R2.fq.gz \
--adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
--adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
--minquality 25 \
--minlength 20 \
--basename "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed \
--trimns \
--trimqualities \
--collapse

# File removal
# rm -f \
# "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.discarded \
# "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.settings \
# "$WD"/"$dataprep"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.singleton.truncated \
# "$WD"/temp/"$(basename "$sample")"_R1.fq.gz \
# "$WD"/temp/"$(basename "$sample")"_R2.fq.gz

exit 0