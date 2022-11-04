#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:10:00

SD="BACKUP/population_genetics/collembola/Orchesella_cincta"

WD="people/Jeppe_Bayer/steps"

# for sample in "$SD"/*; do
    
#     if [ "$(ls -A "$sample")" ]; then
#         for file in "$sample"/*.bam; do
#             if [ ! -e "$file" ]; then
#                 echo "$sample has not been processed"
#                 [[ -d $WD/01_data_preparation/$(basename $SD)/"$(basename "$sample")" ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename $SD)"/"$(basename "$sample")"

#             else
#                 echo "$sample has been processed, $file"
#             fi
#         done
#     else
#         echo "$sample is empty"
#     fi
# done





# [[ -d $WD/01_data_preparation/$(basename $SD) ]] || mkdir -m 764 "$WD"/01_data_preparation/"$(basename $SD)"

# sbatch people/Jeppe_Bayer/scripts/test2.sh $WD





R1=

for R2 in BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F/* ; do
    if [ "$R1" ] ; then
        # 2nd entry - pair
        srun echo "R1 is $R1, R2 is $R2"
    else
        # First entry - just remember.
        R1=$R2
    fi
done


exit 0