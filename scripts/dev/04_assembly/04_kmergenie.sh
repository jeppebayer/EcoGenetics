#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 15G
#SBATCH --cpus-per-task 30
#SBATCH --time 06:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/kmergenie_Enic-%j.out

kmergenie \
/faststorage/project/EcoGenetics/BACKUP/PacBio_HiFi/Entomobrya_nicoleti/ENIC21_m64101e_221211_025112.hifi_reads.fastq.gz \
-l 17 \
-k 30 \
-s 1 \
-o /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/kmergenie_Enic.txt \
-t 30 \
--diploid

exit 0