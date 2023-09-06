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

# Number of bp positions to process at a time
# bpsperround=500
 
nname=""
if [ $((bestn)) -gt 0 ]; then
    nname="_n=$bestn"
fi

# Holds list of all SQ names and LNs from BAM file
SQlist=$(samtools view -H "$targetbam" | rg '@SQ' | awk '{print $2}')
LNlist=$(samtools view -H "$targetbam" | rg '@SQ' | awk '{print $3}')

# Turns SQlist and LNlist into an arrays to easilier access individual SQ names and LNs
x=1
for line in $SQlist; do
    SQarray["$x"]="$line"
    ((x++))
done

x=1
for line in $LNlist; do
    LNarray["$x"]="$line"
    ((x++))
done

# Number of SQ in reference
nSQ=${#SQarray[@]}

# # Number of bp positions to process at a time
# bpsperround=50000

# SQroundsarray=()
# LNroundsarray=()
# numarray=()
# for ((i=1; i<="$nSQ"; i++)); do
#     LN=${LNarray[i]:3}
#     nwhole=$((LN/bpsperround))
#     nremainder=$((LN%bpsperround))
#     if [ "$nwhole" == 0 ]; then
#         SQroundsarray+=(SQlist["$i"])
#         LNroundsarray+=("0-$LN")
#         numarray+=("1")
#     else
#         nrounds=$((nwhole + nremainder))
#         for ((n=1; n<="$nrounds"; n++)); do
#             SQroundsarray+=(SQlist["$i"])
#             numarray+=("$n")
#             if [ "$n" == "$nrounds" ]; then
#                 start=$(((n - 1) * bpsperround))
#                 end=${LN}
#                 LNroundsarray+=("$start-$end")
#             else
#                 start=$(((n - 1) * bpsperround))
#                 end=$((n * bpsperround))
#                 LNroundsarray+=("$start-$end")
#             fi
#         done
#     fi
# done

# arraynum="${SLURM_ARRAY_TASK_ID}"
# region=${SQroundsarray[arraynum]:3}

# echo -e "Current region is $region\nArray number is $arraynum\nSQ is ${SQroundsarray[arraynum]}"
# echo -e "The covered bp postions are ${LNroundsarray[arraynum]}\n"
# # Uses report monomorphic to have all positions, even those with no variation, included. min alternate fraction set to 0 so the only criteria for considering alternate variants is that it is observed at least twice
# freebayes -f "$RG" -n "$bestn" -p "$ploidy" -r "$region":"${LNroundsarray[arraynum]}" --min-alternate-fraction 0 --pooled-discrete --report-monomorphic "$targetbam" > "$temp"/"$sample_name""$nname"_"$region"_"${numarray[arraynum]}".vcf
# # Adds information to out file
# lines=$(rg -c "^[^#]" "$temp"/"$sample_name""$nname"_"$region"_"${numarray[arraynum]}".vcf)
# [ -z "$lines" ] && lines=0
# echo -e "Number of variant lines in VCF is $lines"
# exit 0

##############################################################################################################################

SQnum="${SLURM_ARRAY_TASK_ID}"

# Slices SQ names to work as region name
region=${SQarray[SQnum]:3}

# Creates an array which contains the first and last bp position of the specific region
# IFS=$'\n' read -r -d '' -a bpposition < <(sed -n '1p ; $p' <(samtools view "$targetbam" "$region") | bioawk -c sam '{print $pos}')

# if [ ${#bpposition[@]} == 0 ]; then
#     echo -e "No sample sequence has been mapped to this reference region"
#     exit 0
# fi

# Slices LN to work as bp position name
lastposition=${LNarray[SQnum]:3}
bpposition=(0 "$lastposition")

# Adds information to out file
echo -e "Current region is $region\nSQnum is $SQnum\nSQ is ${SQarray[SQnum]}"

#############################
# Adds information to out file
echo -e "The 'first:last' bp postion is ${bpposition[0]}:${bpposition[1]}\n"
# Uses report monomorphic to have all positions, even those with no variation, included. min alternate fraction set to 0 so the only criteria for considering alternate variants is that it is observed at least twice
freebayes -f "$RG" -n "$bestn" -p "$ploidy" -r "$region" --min-alternate-fraction 0 --min-coverage 200 --pooled-discrete  "$targetbam" > "$temp"/"$sample_name"_"$region""$nname"_complete.vcf
# Adds information to out file
lines=$(rg -c "^[^#]" "$temp"/"$sample_name"_"$region""$nname"_complete.vcf)
[ -z "$lines" ] && lines=0
echo -e "Number of variant lines in VCF is $lines"
exit 0
#############################

# Number of rounds processing n bp positions at a time needed
nrounds=$((bpposition[1]/bpsperround))

# Number of bp positions to run in last round
remainder=$((bpposition[1]%bpsperround))

# If remainder = 1 or 0 then the last round can fail, this avoids that case
if [ $((remainder)) -le 1 ]; then
    nrounds=$((nrounds - 1))
    remainder=$((remainder + bpsperround))
fi

# If region is so short that it wouldn't make for more than 1 whole round whole region is done in one round 
if [ $((nrounds)) -le 1 ]; then
    # Adds information to out file
    echo -e "The 'first:last' bp postion is ${bpposition[0]}:${bpposition[1]}\nVery short region - Will run 1 round with ${bpposition[1]} bp positions\n"
    # Uses report monomorphic to have all positions, even those with no variation, included. min alternate fraction set to 0 so the only criteria for considering alternate variants is that it is observed at least twice
    freebayes -f "$RG" -n "$bestn" -p "$ploidy" -r "$region" --min-alternate-fraction 0 --pooled-discrete --report-monomorphic "$targetbam" > "$temp"/"$sample_name"_"$region""$nname"_complete.vcf
    # Adds information to out file
    lines=$(rg -c "^[^#]" "$temp"/"$sample_name"_"$region""$nname"_complete.vcf)
    [ -z "$lines" ] && lines=0
    echo -e "Number of variant lines in VCF is $lines"
    exit 0
fi

# Adds information to out file
echo -e "The 'first:last' bp postion is ${bpposition[0]}:${bpposition[1]}\nNumber of full rounds with $bpsperround bp positions will be $nrounds\nThe remainder consists of $remainder bp positions for 1 round"

for ((i=0; i<="$nrounds" - 1; i++)); do
    # Creates interval of bp positions for creating VCF file
    start=$(((i * bpsperround) + 1))
    end=$(((bpsperround * i) + bpsperround))

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

    # Adds information to out file
    echo -e "\nStart-end interval for round $num is $start-$end\n"

    # Runs freebayes on subset of bp positions
    # Uses report monomorphic to have all positions, even those with no variation, included. min alternate fraction set to 0 so the only criteria for considering alternate variants is that it is observed at least twice
    freebayes -f "$RG" -n "$bestn" -p "$ploidy" -r "$region":"$start"-"$end" --min-alternate-fraction 0 --pooled-discrete --report-monomorphic "$targetbam" > "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf

    # Adds information to out file
    lines=$(rg -c "^[^#]" "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf)
    [ -z "$lines" ] && lines=0
    echo -e "Number of variant lines in VCF $num is $lines"

    # Concatenates previous and current VCF
    if [ $((i)) -eq 1 ]; then
        bcftools concat -O v -o "$temp"/"$sample_name"_"$region""$nname"_complete_"$num".vcf "$temp"/"$sample_name"_"$region""$nname"_"$prevnum".vcf "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf

        # Removes bcftools concat info from header
        catline=$(rg -n "^##bcftools_concatCommand" "$temp"/"$sample_name"_"$region""$nname"_complete_"$num".vcf | cut -f1 -d:)
        sed -i "${catline}d" "$temp"/"$sample_name"_"$region""$nname"_complete_"$num".vcf

        # Removes separate VCF parts
        rm -f "$temp"/"$sample_name"_"$region""$nname"_"$prevnum".vcf
        rm -f "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf
    elif [ $((i)) -gt 1 ]; then
        bcftools concat -O v -o "$temp"/"$sample_name"_"$region""$nname"_complete_"$num".vcf "$temp"/"$sample_name"_"$region""$nname"_complete_"$prevnum".vcf "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf
        
        # Removes bcftools concat info from header
        catline=$(rg -n "^##bcftools_concatCommand" "$temp"/"$sample_name"_"$region""$nname"_complete_"$num".vcf | cut -f1 -d:)
        sed -i "${catline}d" "$temp"/"$sample_name"_"$region""$nname"_complete_"$num".vcf

        # Removes separate VCF parts
        rm -f "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf
        rm -f "$temp"/"$sample_name"_"$region""$nname"_complete_"$prevnum".vcf
    fi
done

# Creates start of interval containing the remainder of bp positions
start=$((bpposition[1] - remainder + 1))

# Creates number for remainder round
prevnum="$nrounds"
num=$((nrounds + 1))

# Adds information to out file
echo -e "\nStart-end interval for remainder round $num is $start-${bpposition[1]}\n"

# Runs freebayes on the remainder bp positions
# Uses report monomorphic to have all positions, even those with no variation, included. min alternate fraction set to 0 so the only criteria for considering alternate variants is that it is observed at least twice
freebayes -f "$RG" -n "$bestn" -p "$ploidy" -r "$region":"$start"-"${bpposition[1]}" --min-alternate-fraction 0 --pooled-discrete --report-monomorphic "$targetbam" > "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf

# Adds information to out file
lines=$(rg -c "^[^#]" "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf)
[ -z "$lines" ] && lines=0
echo -e "Number of variant lines in VCF $num is $lines"

# Concatenates previous rounds and remainder
bcftools concat -O v -o "$temp"/"$sample_name"_"$region""$nname"_complete.vcf "$temp"/"$sample_name"_"$region""$nname"_complete_"$prevnum".vcf "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf

# Removes bcftools concat info from header
catline=$(rg -n "^##bcftools_concatCommand" "$temp"/"$sample_name"_"$region""$nname"_complete.vcf | cut -f1 -d:)
sed -i "${catline}d" "$temp"/"$sample_name"_"$region""$nname"_complete.vcf

# Removes separate VCF parts
rm -f "$temp"/"$sample_name"_"$region""$nname"_complete_"$prevnum".vcf
rm -f "$temp"/"$sample_name"_"$region""$nname"_"$num".vcf

exit 0