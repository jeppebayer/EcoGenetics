#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 15G
#SBATCH --cpus-per-task 32
#SBATCH --time 06:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/kmer_genome_est-%j.out

jellyfish count \
-t 32 \
-C \
-m \ # length of kmer
-s 8G \
-o /faststorage/peoject/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/jellyfish_Enic \
/faststorage/project/EcoGenetics/BACKUP/PacBio_HiFi/Entomobrya_nicoleti/ENIC21_m64101e_221211_025112.hifi_reads.fastq.gz