#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 15G
#SBATCH --cpus-per-task 32
#SBATCH --time 10:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/Lela_hifiasm-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

hifiasm \
-o /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/Lepidocyrtus_sp/Lela.asm \
-t 32 \
-k 27
/faststorage/project/EcoGenetics/BACKUP/PacBio_HiFi/Lepidocyrtus_sp/LELA_m64101e_221212_121647.hifi_reads.fastq.gz

exit 0