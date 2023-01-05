#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00

# # Directory containing scripts (Do NOT end with '/')
# scripts="people/Jeppe_Bayer/scripts"

# # Species specific reference genome (in FASTA format)
# RG="BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna"

# # Species specific sample directory (Do NOT end with '/')
# SD="BACKUP/population_genetics/collembola/Orchesella_cincta"

# # Working directory (Do NOT end with '/')
# WD="people/Jeppe_Bayer/steps"

# usage()
# {
# cat << EOF

# Usage: 02_00_init_data_prep.sh [-r|--reference] <reference_genome> [-s|--species] <species_directory> [-d|--directory] <working_directory> [-a|--algorithm] <algorithm> [-h|--help]

# This script is used for initializing the standardized data preparation procedure for sequence data

# PARAMETERS:
#     -r | --reference    Species specific reference genome, abosolute path (reference genome in FASTA format)
#     -s | --species      Species specific sample directory, abosolute path (Do NOT end with '/')
#     -d | --directory    Working directory, abosolute path (Do NOT end with '/')

# OPTIONS:
#     -a | --algorithm    Choice of algorithm to be used during alignment. mem (>70MB, contemporary samples)[default] or aln (<70MB, historic samples)
#     -h | --help         Show this message

# EOF
# }

# RG="/home/jepe/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1/GCA_001718145.1_ASM171814v1_genomic.fna"

# for file in "$(dirname "$RG")"/*.gff; do
#     gff=$file
# done

# echo "$gff"

# Gets path to script location
# If run through sbatch or srun:
# if [ -n "$SLURM_JOB_ID" ]; then
#     script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
#     script_path=$(dirname "$script_path")
# # If run on the frontend:
# else
#     script_path=$(realpath "$0")
# fi

# echo "$script_path"

# count=$(find $SD/Ocin_NYS-F -maxdepth 1 -type f -name '*.fq.gz' | wc -l)
# if [ "$count" == 1 ] ; then
#   echo "$count is 1"
# else
#   echo "$count is not 1"
# fi

# exit 1

# path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# path=$(dirname "$path")

# touch "$path"/testlog.txt

# # Algorithm to use during alignment: mem (>70MB, contemporary samples) [default] or aln (<70MB, historic samples)
# algo="mem"

# # Define memory per cpu in G (must be integer)
# memory="8"

# # Define number of cpus to be used (must be integer)
# cpus="8"

# all_variables=("$RG" "$SD" "$WD" "$algo" "$memory" "$cpus")

# for ref in /home/"$USER"/EcoGenetics/BACKUP/reference_genomes/"$(basename "$SD")"/*.fna; do
#     RG=$ref
#     break
# done

# AdapterRemoval
# jid1=$(sbatch --parsable "$path"/test2.sh > "$path"/testlog.txt)

# Aligning to reference
# sbatch /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/test3.sh "${all_variables[@]}"
# sbatch /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/test3.sh "$RG"

# queue()
# {
# cat << EOF
# #!/bin/bash
# #SBATCH --account EcoGenetics
# #SBATCH --partition normal
# #SBATCH --mem-per-cpu $1G
# #SBATCH --cpus-per-task $2
# #SBATCH --time 00:30:00
# #SBATCH --output=isit3.out

# echo "Hello"

# exit 0
# EOF
# }

# SD="this/path/is/no/good/"

# if [ "${SD: -1}" == "/" ]; then
#     SD=${SD:0:((${#SD} - 1))}
# fi

# echo $SD

# lines=$(wc -l <people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_qname.txt)
# echo "$lines"
# filesize=$(wc -c <BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F/Ocin_NYS-F.pileup)
# echo "$filesize"
# sizeratio=$(awk -v filesize="$filesize" 'BEGIN { print ( filesize / 234060585564) }')
# echo "$sizeratio"
# lineratio=$(awk -v lines="$lines" 'BEGIN { print ( 9402 / lines ) }')
# echo "$lineratio"
# adjustment=$(awk -v sizeratio="$sizeratio" -v lineratio="$lineratio" 'BEGIN { print ( sizeratio * lineratio ) }')
# echo "$adjustment"
# adjbase=$(awk -v adjustment="$adjustment" 'BEGIN { print int( 600 * adjustment ) }')
# echo "$adjbase"
# if [ "$adjbase" -lt 240 ]; then 
#     awk -v adjbase="$adjbase" 'BEGIN { print int( adjbase + 240 ) }'
# else
#     echo "$adjbase"
# fi

# target="BACKUP/PacBio_HiFi/Entomobrya_nicoleti/ENIC21_m64101e_221211_025112.hifi_reads.fastq.gz"
# name=$(basename "$target")
# firstletter=${name:0:1}
# nextletters="${name:1:3}"
# lowerletters="${nextletters,,}"
# echo "$firstletter$lowerletters"

WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Isotoma_sp/purge_dups/01"

# Creates files with all sequence lengths sorted longest to shortest
bioawk \
-v OFS='\t' \
-c fastx \
'{ print $name, length($seq) }' \
"$WD"/purged.fa \
| sort -k1,1 -k2,2nr \
> "$WD"/sequence_lengths

# Calculates total assembly length
total=$(awk '{Total=Total+$2} END{print Total}' "$WD"/sequence_lengths)
awk '{Total=Total+$2} END{print "Total assembly length is: " Total}' "$WD"/sequence_lengths > "$WD"/asm_lenght_"$total"

exit 0