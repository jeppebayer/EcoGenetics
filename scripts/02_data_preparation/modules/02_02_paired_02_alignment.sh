#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Concatenating collapsed single-end files
cat \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.collapsed \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.collapsed.truncated \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.all_collapsed

# Align sample to reference genome
bwa "$algo" -t "$cpus" \
"${RG%.*}" \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.all_collapsed \
| \

# Sort with regards to QNAME and convert to bam format
samtools sort -@ "$(("$cpus" - 1))" -n -O BAM \
-T "$WD"/temp/  \
-o "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_collapsed_aligned.bam \
-

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.collapsed \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.collapsed.truncated \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.all_collapsed

exit 0