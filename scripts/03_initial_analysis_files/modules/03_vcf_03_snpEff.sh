#!/bin/bash

WD="$1" # Working directory
sample_name="$2" # Name of current sample
bestn="$3" # Best n alleles value
dbgenome="$4" # Name of genome in snpEff database

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

nname=""
if [ $((bestn)) -gt 0 ]; then
    nname="_n=$bestn"
fi

sample_vcf="$WD"/"$sample_name""$nname"_complete.vcf
sample_vcf_noext="$(basename "${sample_vcf%.*}")"

# snpEff 
snpeff_path=("$(dirname "$(dirname "$(which python)")")"/share/snpeff*)
snpeffconfig="${snpeff_path[*]}"/snpEff.config

snpEff ann -csvStats "$WD"/"$sample_name"_snpEff_summary.csv -c "$snpeffconfig" -v "$dbgenome" "$sample_vcf" > "$WD"/"$sample_vcf_noext".ann.vcf

exit 0