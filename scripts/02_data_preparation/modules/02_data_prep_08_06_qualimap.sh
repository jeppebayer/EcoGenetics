#!/bin/bash

RG="$1" # Reference genome
SD="$2" # Species directory
sample="$3" # Sample directory

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate data_prep
fi

shopt -s extglob

# Change of JAVA_OPTS in qualimap script
qualimap_path=$(dirname "$(which python)")/qualimap
grep -qxF '\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE"' "$qualimap_path" || sed -i '47s#.*#\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE"#' "$qualimap_path"

# First checks whether a .gff file for the reference genome is available
for file_gff in "$(dirname "$RG")"/*.gff; do
    if [ -e "$file_gff" ]; then
        # .gff file is available
        gff="$file_gff"

        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -gff "$gff" \
        -outdir "$SD"/"$(basename "$sample")"/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=20G
        exit 0
    else
        for file_gtf in "$(dirname "$RG")"/*.gtf; do
            if [ -e "$file_gtf" ]; then
                # .gtf file is available
                gtf=$file_gtf

                qualimap bamqc \
                -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
                -gff "$gtf" \
                -outdir "$SD"/"$(basename "$sample")"/qualimap \
                -outfile "$(basename "$sample")"_qualimap.pdf \
                -outformat PDF \
                --java-mem-size=20G
                exit 0
            else
                # Neither .gff nor .gtf file is not available
                qualimap bamqc \
                -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
                -outdir "$SD"/"$(basename "$sample")"/qualimap \
                -outfile "$(basename "$sample")"_qualimap.pdf \
                -outformat PDF \
                --java-mem-size=20G
                exit 0
            fi
        done
    fi
done
