#!/bin/bash

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

dir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/03_initial_analysis_files/Entomobrya_nicoleti

out="$dir"/out
arraynums=""

for folder in "$dir"/*; do

    for file in "$folder"/out/*; do

        if rg -Fq "CANCELLED AT" "$file"; then
            name="$(basename "$file")"
            name=${name%-*.out}
            name=${name#*vcf-*-*}
            arraynums+="$name",
            # rm "$file"
        fi
    done
    echo "$(basename "$folder")"
    echo "$arraynums"
done

exit 0