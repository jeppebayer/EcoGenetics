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

# Change of JAVA_OPTS in qualimap script
qualimap_path=$(dirname "$(which python)")/qualimap
grep -qxF '\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE"' "$qualimap_path" || sed -i '47s#.*#\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE"#' "$qualimap_path"

# First checks whether a .gff file for the reference genome is available
for file in "$(dirname "$RG")"/*.gff; do
    if [ -e "$file" ]; then
        
        # .gff file is available
        gff=$file

        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -gff "$gff" \
        -outdir "$SD"/"$(basename "$sample")"/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=20G
        exit 0

    else

        # .gff file is not available
        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -outdir "$SD"/"$(basename "$sample")"/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=20G
        exit 0
    fi
done
