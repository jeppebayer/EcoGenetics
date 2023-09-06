#!/bin/bash

reference_genome="$1"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate repeatmasking
fi

BuildDatabase \
-name "$(basename "$(dirname "$reference_genome")")" \
-engine rmblast \
"$reference_genome"

exit 0