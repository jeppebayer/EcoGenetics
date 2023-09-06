#!/bin/bash

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

if [ -z "$1" ]; then
    configuration_name=$(readlink -f "template")
else
    configuration_name="$(readlink -f "$1")"
fi

configuration_content() {
cat << EOF
#Configuration file for:
#03_initial_analysis_files  :   $(basename "$configuration_name")

#Paths to files and folders should be absolute.
#Comments can be added by starting the line with '#'. Starting lines with '>' is reserved for settings.

>Species name:
[genus] [species]

>Reference genome:
[path/to/genome.fasta]

#One or more BAM files containing sample alignments to be analyzed.
#If using multiple files add one file per line.
>BAM:
[path/to/sample.bam]

#In the given working directory species and sample specific folders will be made automatically.
>Working directory:
$(dirname "$configuration_name")

#Ploidy of the sample(s).
#When samples are pooled it should equal the number of genome copies in the pool.
>Ploidy:
100

#Freebayes parameter 'use-best-n-alleles'.
#Sets the number of haplotypes to consider for a given position.
>use-best-n-alleles':
3

#The database entry can be chosen from the custom data entries using any entry under 'Genome_name'.
#The custom database entries can be view at the bottom of the configuration file.
#Any standard entry in Snpeff can also be used.
>Database Entry
[Genome_name]

#Job steps.
#The presence or absence of an 'x' below each numbered step determines whether or not that step is to be executed.
>1. Create VCF parts
x

>2. Concatenate all VCF parts
x

>3. Annotation of VCF using snpEff
x

>4. Create SFS
x

#-----CUSTOM DATABASE ENTRIES-----
#As of $(date +"%A%t%H:%M%t%d-%m-%Y")

$db_status
EOF
}

snpeff_path=("$(dirname "$(dirname "$(which python)")")"/share/snpeff*)
custom_dbs="${snpeff_path[*]}"/customDBs.txt
if [ -e "$custom_dbs" ]; then
    db_status="$(awk -F "\t" '{printf "%-35s %s\n", "#"$1, $2}' "$custom_dbs")"
elif [ ! -e "${snpeff_path[*]}" ]; then
    db_status="#snpEff does not appear to be installed in the current environment"
else
    db_status="#There seems to be no custom snpEff database entries"
fi

configuration_content > "$configuration_name".03iaf.configuration

exit 0