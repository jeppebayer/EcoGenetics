#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 02:00:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/hifiadapterfilt-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly

fi

export PATH=$PATH:/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/HiFiAdapterFilt
export PATH=$PATH:/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/HiFiAdapterFilt/DB

# in_file="$(readlink -f "$1")"
# WD="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly"

# bash /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/HiFiAdapterFilt/pbadapterfilt.sh \
# -p "$(basename "${in_file%.*.*}")" \
# -t 10 \
# -o "$WD"/"$(basename "$(dirname "$in_file")")"

bash /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/HiFiAdapterFilt/pbadapterfilt.sh \
-t 10

exit 0