#!/bin/bash

targetbam="$1" # Target bam file
RG="$2" # Reference genome
temp="$3" # Temp directory
sample_name="$4" # Name of current sample
bestn="$5" # Best n alleles value
ploidy="$6" # Ploidy value for sample

export TMPDIR="$temp"

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

repeat_bed=("$(dirname "$RG")"/*.repeatregions.bed)

nname=""
if [ $((bestn)) -gt 0 ]; then
    nname="_n=$bestn"
fi

# Holds list of all SQ names and LNs from BAM file
SQlist=$(samtools view -H "$targetbam" | rg '@SQ' | awk '{print $2}')
LNlist=$(samtools view -H "$targetbam" | rg '@SQ' | awk '{print $3}')

# Slices and turns SQlist and LNlist into arrays to easilier access individual SQ names and LNs
x=1
for line in $SQlist; do
    SQarray["$x"]="${line:3}"
    ((x++))
done

x=1
for line in $LNlist; do
    LNarray["$x"]="${line:3}"
    ((x++))
done

SQnum="${SLURM_ARRAY_TASK_ID}"

# Access region name
region=${SQarray[SQnum]}

# Access last bp position name
lastposition=${LNarray[SQnum]}

# Adds information to out file
echo -e "Current region is $region\nSQnum is $SQnum\nSQ is ${SQarray[SQnum]}"
echo -e "The 'first:last' bp position is 0:$lastposition\n"

# Min alternate fraction set to 0 so the only criteria for considering alternate variants is that it is observed at least twice
# freebayes -f "$RG" -n "$bestn" -p "$ploidy" -r "$region" --min-alternate-fraction 0 --min-coverage 200 --pooled-discrete  "$targetbam" > "$temp"/"$sample_name"_"$region""$nname"_complete.vcf

# Current version run on OrcCin-NYS-F
freebayes \
-f "$RG" \
-n "$bestn" \
-p "$ploidy" \
-r "$region" \
--min-alternate-fraction 0 \
--min-coverage 200 \
--pooled-discrete  \
"$targetbam" \
> "$temp"/"$sample_name"_"$region""$nname"_complete.vcf

# # Version with filtering build in 
# freebayes \
# -f "$RG" \
# -n "$bestn" \
# -p "$ploidy" \
# -r "$region" \
# --min-alternate-fraction 0 \
# --report-monomorphic \
# --min-coverage 200 \
# --pooled-discrete  \
# "$targetbam" \
# | SnpSift intervals \
# -x \
# "${repeat_bed[*]}" \
# | SnpSift filter \
# "((QUAL >= 30) | (QUAL = 0)) && (TYPE = 'snp')" \
# > "$temp"/"$sample_name"_"$region""$nname"_complete.vcf

# Adds information to out file
lines=$(rg -c "^[^#]" "$temp"/"$sample_name"_"$region""$nname"_complete.vcf)
[ -z "$lines" ] && lines=0
echo -e "Number of variant lines in VCF is $lines"

exit 0