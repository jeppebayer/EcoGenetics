#!/bin/bash

# # Rename sample directories
# for dir in /faststorage/project/EcoGenetics/BACKUP/population_genetics/*/*/*; do
#     # Avoids GERP directory
#     if [[ "$dir" == */GERP/* ]]; then
#         continue
#     fi
#     # Creates species abbreviation
#     species_name="$(basename "$(dirname "$dir")")"
#     genus=${species_name%_*}
#     genus=${genus::3}
#     genus=${genus^}
#     species=${species_name#*_}
#     species=${species::3}
#     species=${species^}
#     abbr="$genus""$species"

#     # Gets location information of original name
#     current_dirname="$(basename "$dir")"
#     location=${current_dirname#*_}

#     # Renames directory
#     mv "$dir" "$(dirname "$dir")"/"$abbr"_"$location"
# done

# # Rename sample files
# for file in /faststorage/project/EcoGenetics/BACKUP/population_genetics/*/*/*/*; do
#     # Avoids GERP directory
#     if [[ "$file" == */GERP/* ]]; then
#         continue
#     fi
#     if [[ "$file" == *.bam ]] || [[ "$file" == *.bam.bai ]]; then
#         # Creates species abbreviation
#         species_name="$(basename "$(dirname "$(dirname "$file")")")"
#         genus=${species_name%_*}
#         genus=${genus::3}
#         genus=${genus^}
#         species=${species_name#*_}
#         species=${species::3}
#         species=${species^}
#         abbr="$genus""$species"

#         # Gets location information of original name
#         current_filename="$(basename "$file")"
#         namepart=${current_filename#*_}

#         # Renames file
#         mv "$file" "$(dirname "$file")"/"$abbr"_"$namepart"
#     fi
# done

# # Rename Qualimap files
# for file in /faststorage/project/EcoGenetics/BACKUP/population_genetics/*/*/*/qualimap/*; do
#     # Avoids GERP directory
#     if [[ "$file" == */GERP/* ]]; then
#         continue
#     fi
#     if [[ "$file" == *.pdf ]]; then
#         # Creates species abbreviation
#         species_name="$(basename "$(dirname "$(dirname "$(dirname "$file")")")")"
#         genus=${species_name%_*}
#         genus=${genus::3}
#         genus=${genus^}
#         species=${species_name#*_}
#         species=${species::3}
#         species=${species^}
#         abbr="$genus""$species"

#         # Gets location information of original name
#         current_filename="$(basename "$file")"
#         namepart=${current_filename#*_}

#         # Renames file
#         mv "$file" "$(dirname "$file")"/"$abbr"_"$namepart"
#     fi
# done

data_prep=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/02_data_preparation

for speciesdir in "$data_prep"/*; do
    # Creates species abbreviation from species name
    species_name=$(basename "$speciesdir")
    genus=${species_name%_*}; genus=${genus::3}; genus=${genus^}
    species=${species_name#*_}; species=${species::3}; species=${species^}
    speciesabbr="$genus""$species"
    for sampledir in "$speciesdir"/*; do
        samplename=$(basename "$sampledir")
        if [ "$samplename" == "EntNic_JEJ-C119" ] || [ "$samplename" == "EntNic_ULJ-C122" ]; then
            continue
        fi
        samplenamepart=${samplename#*_}
        mv "$sampledir" "$(dirname "$sampledir")"/"$speciesabbr"_"$samplenamepart"
        for bam in "$sampledir"/*; do
            if [ ! -d "$bam" ]; then
                bamname=$(basename "$bam")
                bamnamepart=${bamname#*_}
                mv "$bam" "$(dirname "$bam")"/"$speciesabbr"_"$bamnamepart"
            elif [ -d "$bam" ]; then
                for file in "$bam"/*; do
                    filename=$(basename "$file")
                    filenamepart=${filename#*_}
                    mv "$file" "$(dirname "$file")"/"$speciesabbr"_"$filenamepart"
                done
            fi
        done
    done
done

exit 0