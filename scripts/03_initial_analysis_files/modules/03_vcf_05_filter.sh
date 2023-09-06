#!/bin/bash

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

# QUALITY >= 30 (accuracy 99.9%)
# Coverage >= 25

sample_vcf=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_initial_analysis_files/Orchesella_cincta/OrcCin_NYS-F/OrcCin_NYS-F_n=3_complete.ann.vcf

SnpSift filter -f $sample_vcf "(QUAL >= 30) && (DP >= 25) && (TYPE = 'snp') && ((ANN[*].EFFECT has 'intergenic_region') | (ANN[*].EFFECT has 'conserved_intergenic_variant'))" > "$(dirname "$sample_vcf")"/filtered.vcf

exit 0