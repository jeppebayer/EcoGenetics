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

# Add fixmate tag to alignment. Can only be done on name sorted alignment
samtools fixmate -@ "$(("$cpus" - 1))" -m -O BAM \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
- | \

# Position sort alignment and pipe output
samtools sort -@ "$(("$cpus" - 1))" -O BAM \
-T "$WD"/temp/ \
- | \

# Mark duplicates and save output as .coordinatesort.bam. Also outputs some realted statistics
samtools markdup -@ "$(("$cpus" - 1))" -s \
-f "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.markdupstats \
-T "$WD"/temp/ \
- \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam

exit 0