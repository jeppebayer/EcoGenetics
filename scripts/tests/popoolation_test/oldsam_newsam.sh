#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem 12G
#SBATCH --cpus-per-task 2
#SBATCH --time 04:00:00
#SBATCH --chdir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/popoolation_test
#SBATCH --output=pileup_test-%j-%A-%a.out
#SBATCH --array=1-3

input_file=/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/OrcCin_NYS-F/OrcCin_NYS-F_filtered.bam
reference_genome=/faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna

if [ "${SLURM_ARRAY_TASK_ID}" == "1" ]; then

    # Sources necessary environment
    if [ "$USER" == "jepe" ]; then
        # shellcheck disable=1090
        source /home/"$USER"/.bashrc
        #shellcheck disable=1091
        source activate oldsam
    fi

    echo "I'm Old SAMtools, aka classic pileup"

    samtools pileup \
        -f "$reference_genome" \
        "$input_file" \
        > OrcCin-NYS-F.oldsam.pileup

elif [ "${SLURM_ARRAY_TASK_ID}" == "2" ]; then

    # Sources necessary environment
    if [ "$USER" == "jepe" ]; then
        # shellcheck disable=1090
        source /home/"$USER"/.bashrc
        #shellcheck disable=1091
        source activate data_prep
    fi

    echo "I'm New SAMtools, aka mpileup"

    samtools mpileup \
        -f "$reference_genome" \
        -o  OrcCin-NYS-F.newsam.pileup \
        "$input_file"

elif [ "${SLURM_ARRAY_TASK_ID}" == "3" ]; then

    # Sources necessary environment
    if [ "$USER" == "jepe" ]; then
        # shellcheck disable=1090
        source /home/"$USER"/.bashrc
        #shellcheck disable=1091
        source activate vcf
    fi
    
    echo "I'm BCFtools, aka what happened to classic pileup"

    bcftools mpileup \
        -f "$reference_genome" \
        -o OrcCin-NYS-F.bcftools.pileup \
        "$input_file"

fi

exit 0