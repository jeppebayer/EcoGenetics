#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:10:00

# Indexing reference
# jid1=$(sbatch --parsable /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/01.sh)

# Align sample to reference genome
### Use if reference genome is NOT indexed ###
# jid2=$(sbatch --parsable --dependency=afterany:"$jid1" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/02.sh)
### Use if reference genome IS indexed ###
# jid2=$(sbatch --parsable /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/02.sh)

# Convert sam format to bam format
# jid3=$(sbatch --parsable --dependency=afterany:"$jid2" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/03.sh)

# Sort alignment with regard to QNAME
# jid4=$(sbatch --parsable --dependency=afterany:"$jid3" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/04.sh)

# Add fixmate tag to alignment
# jid5=$(sbatch --parsable --dependency=afterany:"$jid4" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/05.sh)

# jid5=$(sbatch --parsable /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/05.sh)


# Position sort alignment
# jid6=$(sbatch --parsable --dependency=afterany:"$jid5" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/06.sh)

# Mark duplicates
# jid7=$(sbatch --parsable --dependency=afterany:"$jid6" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/07.sh)

# jid7=$(sbatch --parsable /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/07.sh)

# Creates bai index for alignment
# jid8=$(sbatch --parsable --dependency=afterany:"$jid7" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/08.sh)

# Creates flagstat file for alignment 
# jid9=$(sbatch --parsable --dependency=afterany:"$jid8" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/09.sh)

# Creates idxstats file for alignment
# jid10=$(sbatch --parsable --dependency=afterany:"$jid9" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/10.sh)

# Creates coverage file for alignment
# jid11=$(sbatch --parsable --dependency=afterany:"$jid10" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/11.sh)

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
# jid12=$(sbatch --parsable --dependency=afterany:"$jid11" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/12.sh)

jid12=$(sbatch --parsable /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/12.sh)

# Creates files listing number of reads in markdup file on first line
# number of reads in filtered file on second line
# and % of remaining reads on the third line
sbatch --dependency=afterany:"$jid12" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/13.sh

exit 0