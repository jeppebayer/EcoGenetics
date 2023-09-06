#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem 12G
#SBATCH --cpus-per-task 2
#SBATCH --time 08:00:00
#SBATCH --chdir=/faststorage/project/EcoGenetics
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/popoolation_test/pileup-%A_%a.out
#SBATCH --array=1-3

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate oldsam
fi

if [ "${SLURM_ARRAY_TASK_ID}" == "1" ]; then

    input_file=/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/OrcCin_NYS-F/OrcCin_NYS-F_filtered.bam
    reference_genome=/faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna
    output_file=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/popoolation_test/OrcCin_NYS-F.pileup

elif [ "${SLURM_ARRAY_TASK_ID}" == "2" ]; then

    input_file=/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/OrcCin_RYJ/OrcCin_RYJ_filtered.bam
    reference_genome=/faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna
    output_file=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/popoolation_test/OrcCin_RYJ.pileup

elif [ "${SLURM_ARRAY_TASK_ID}" == "3" ]; then

    input_file=/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/OrcCin_NYS-F-C25/OrcCin_NYS-F-C25_filtered.bam
    reference_genome=/faststorage/project/EcoGenetics/BACKUP/reference_genomes/Orchesella_cincta/GCA_001718145.1_ASM171814v1_genomic.fna
    output_file=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/popoolation_test/OrcCin_NYS-F-C25.pileup

fi

samtools pileup \
    -f "$reference_genome" \
    "$input_file" \
    > "$output_file"

exit 0