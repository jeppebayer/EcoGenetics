#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 15G
#SBATCH --cpus-per-task 32
#SBATCH --time 10:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/Ecin_hifiasm-%j.out

hifiasm \
-o Ecin.asm \
-t 32 \
/faststorage/project/EcoGenetics/BACKUP/PacBio_HiFi/Entomobrya_nicoleti/ENIC21_m64101e_221211_025112.hifi_reads.fastq.gz