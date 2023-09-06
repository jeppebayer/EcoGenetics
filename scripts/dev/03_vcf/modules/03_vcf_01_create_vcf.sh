#!/bin/bash

targetbam="$1" # Target bam file
RG="$2" # Reference genome
temp="$3" # Temp directory
sample_name="$4" # Name of current sample

export TMPDIR="$temp"

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

# Holds list of all contig names from bam file
SQlist=$(samtools view -H "$targetbam" | rg '@SQ' | awk '{print $2}')

# Turns SQlist into an array to easilier access individual contig names
x=1
for line in $SQlist; do
    SQarray["$x"]="$line"
    ((x++))
done

SQnum="${SLURM_ARRAY_TASK_ID}"

# Slices contig names to work as region name
region=${SQarray[SQnum]:3}

# Number of records within a given contig
nrecords=$(samtools view -c "$targetbam" "$region")

# Number of rounds processing 5500 records at a time needed
nrounds=$((nrecords/5000))

# Number of records to run in last round
remainder=$((nrecords%5000))

for ((i=0; i<="$nrounds" - 1; i++)); do
    # Creates a subset of records for creating VCF file
    start=$((i * 5000 + 1))
    end=$((5000 * i + 5000))
    recordset=$(sed -n "'$start','$end'p" <(samtools view "$targetbam" "$region"))

    # Adjusts naming according to file number
    round=$((i + 1))
    roundlength=${#round}
    maxlength=${#nrounds}

    if [ $((roundlength)) -lt $((maxlength)) ]; then
        dif=$((maxlength - roundlength))
        prevnum="$num"
        num=""
        for (( j=1;  j<="$dif"; j++ )); do
            num="0$num"
        done
        num="$num$round"
    else
        prevnum="$num"
        num="$round"
    fi

    # Runs freebayes on set of records
    freebayes -f "$RG" -p 100 --pooled-discrete <(samtools view -O BAM <(cat <(samtools view -H "$targetbam") <(echo "$recordset"))) > "$temp"/"$sample_name"_"$region"_"$num".vcf

    # Concatenates previous and current VCF
    if [ $((i)) -eq 1 ]; then
        bcftools concat -O v -o "$temp"/"$sample_name"_"$region"_complete.vcf "$temp"/"$sample_name"_"$region"_"$prevnum".vcf "$temp"/"$sample_name"_"$region"_"$num".vcf
        
        # Removes separate VCF parts
        rm -f "$temp"/"$sample_name"_"$region"_"$prevnum".vcf
        rm -f "$temp"/"$sample_name"_"$region"_"$num".vcf
    elif [ $((i)) -gt 1 ]; then
        bcftools concat -O v -o "$temp"/"$sample_name"_"$region"_complete.vcf "$temp"/"$sample_name"_"$region"_complete.vcf "$temp"/"$sample_name"_"$region"_"$num".vcf
        
        # Removes separate VCF parts
        rm -f "$temp"/"$sample_name"_"$region"_"$num".vcf
    fi
done

# Creates subset containing the remainder of records
start=$((nrecords - remainder + 1))
remainderset=$(sed -n "'$start','$nrecords'p" <(samtools view "$targetbam" "$region"))

# Creates number for remainder round
num=$((nrounds + 1))

# Runs freebayes on the set of remainders
freebayes -f "$RG" -p 100 --pooled-discrete <(samtools view -O BAM <(cat <(samtools view -H "$targetbam") <(echo "$remainderset"))) > "$temp"/"$sample_name"_"$region"_"$num".vcf

# Concatenates previous rounds and remainder
bcftools concat -O v -o "$temp"/"$sample_name"_"$region"_complete.vcf "$temp"/"$sample_name"_"$region"_complete.vcf "$temp"/"$sample_name"_"$region"_"$num".vcf

# Removes separate VCF parts
rm -f "$temp"/"$sample_name"_"$region"_"$num".vcf

exit 0