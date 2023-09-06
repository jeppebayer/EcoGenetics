#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:10:00

# Indexing reference
# jid1=$(sbatch --parsable /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/01.sh)

### Use if reference genome is NOT indexed ###
# jid2=$(sbatch --parsable --dependency=afterany:"$jid1" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/02.sh)
### Use if reference genome IS indexed ###
jid2=$(sbatch --parsable /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/02_alt.sh)

jid3=$(sbatch --parsable --dependency=afterany:"$jid2" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/03_alt.sh)

sbatch --dependency=afterany:"$jid3" /home/jepe/EcoGenetics/people/Jeppe_Bayer/scripts/data_preparation/04_alt.sh

exit 0