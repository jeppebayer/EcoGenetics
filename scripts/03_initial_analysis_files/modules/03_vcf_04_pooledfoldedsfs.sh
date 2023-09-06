#!/bin/bash

WD="$1" # Working directory
sample_name="$2" # Name of current sample
bestn="$3" # Best n alleles value
ploidy="$4" # Ploidity of sample

foldedploidy=$((ploidy / 2))

nname=""
if [ $((bestn)) -gt 0 ]; then
    nname="_n=$bestn"
fi

sample_vcf="$WD"/"$sample_name""$nname"_complete.ann.vcf

if [ ! -e "$sample_vcf" ]; then
    sample_vcf="$WD"/"$sample_name""$nname"_complete.vcf
    if [ ! -e "$sample_vcf" ]; then
        echo "Can't locate VCF file in $WD"
        exit 1
    fi
fi

# First, removes all lines staring with '#'
# Second, removes everything that comes after ':'
# Third, removes all lines containing numbers other than '0' or '1'
# Fourth, changes the field seperator from '/' to '\t'
# Fifth, sum all numbers in a row, then adds 1 to the corresponding index in 'sfsarray', 
# if sum is greater than 50 the index becomes '100 - sum'. Lastly prints all values in 'sfsarray' 
# as a single line with tab separated values
# awk -F "\t" -v foldedploidy=$foldedploidy '
# BEGIN {
#     for (i=0; i<=foldedploidy; i++) {
#         sfsarray[i] = 0
#     }
# }
# {
# sum = 0 ;
# for (i=1; i<=NF; i++) {
#     sum += $i
#     } ;
# if (sum > foldedploidy) {
#     foldsum = (100 - sum)
#     sfsarray[foldsum] += 1
#     } ;
# if (sum <= foldedploidy) {
#     sfsarray[sum] += 1
#     }
# }
# END {
#     sep = "\t"
#     sfs = sfsarray[0]
#     for (n=1; n<=foldedploidy; n++) {
#         sfs = sfs sep sfsarray[n]
#     } ;
#     print sfs
# }' \
#     <(awk -F "/" '
#         BEGIN {
#             OFS = "\t"
#         } 
#         {$1=$1 ; print $0}' \
#             <(awk -F " " '
#                 $0 !~ /[23456789]/ {
#                     print $0
#                     }' \
#                 <(awk -F ":" '
#                 {print $1}' \
#                     <(awk -F "\t" '
#                     /^[^#]/ {
#                         print $NF
#                         }' "$sample_vcf")))) \
# > "$WD"/"$sample_name"_pooled_folded.sfs

awk \
    -F "\t" \
    '/^[^#]/ {
        print $NF
    }' \
    "$sample_vcf" \
| awk \
    -F ":" \
    '{print $1}' \
| awk \
    -F " " \
    '$0 !~ /[23456789]/ {
        print $0
    }' \
| awk \
    -F "/" \
    'BEGIN {
        OFS = "\t"
    } 
    {$1=$1 ; print $0}' \
| awk \
    -F "\t" \
    -v foldedploidy=$foldedploidy \
    'BEGIN {
        for (i=0; i<=foldedploidy; i++) {
            sfsarray[i] = 0
        }
    }
    { sum = 0 ;
    for (i=1; i<=NF; i++) {
        sum += $i
    } ;
    if (sum > foldedploidy) {
        foldsum = (100 - sum)
        sfsarray[foldsum] += 1
    } ;
    if (sum <= foldedploidy) {
        sfsarray[sum] += 1
    } }
    END {
        sep = "\t"
        sfs = sfsarray[0]
        for (n=1; n<=foldedploidy; n++) {
            sfs = sfs sep sfsarray[n]
        } ;
        print sfs
    }' \
    > "$WD"/"$sample_name"_pooled_folded.sfs

exit 0