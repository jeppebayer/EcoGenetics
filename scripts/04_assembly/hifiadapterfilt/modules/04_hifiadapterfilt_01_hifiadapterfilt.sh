#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

targetprefix="$1"
filtdir="$2"
script_path="$3"
currentuser="$4"
originalprefix="$(basename "$5")"

if [ "$USER" == "$currentuser" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate genome_assembly
fi

export PATH=$PATH:"$script_path"/modules/HiFiAdapterFilt
export PATH=$PATH:"$script_path"/modules/HiFiAdapterFilt/DB

bash "$script_path"/modules/HiFiAdapterFilt/pbadapterfilt.sh \
-t 10 \
-p "$(basename "$targetprefix")" \
-o "$(basename "$filtdir")"

if [ "$(basename "$targetprefix")" == "temp_name" ]; then
    mv -f "$filtdir"/"$(basename "$targetprefix")".contaminant.blastout "$filtdir"/"$originalprefix".contaminant.blastout
    mv -f "$filtdir"/"$(basename "$targetprefix")".blocklist "$filtdir"/"$originalprefix".blocklist
    mv -f "$filtdir"/"$(basename "$targetprefix")".filt.fastq.gz "$filtdir"/"$originalprefix".filt.fastq.gz
    mv -f "$filtdir"/"$(basename "$targetprefix")".stats "$filtdir"/"$originalprefix".stats
    rm -f "$(dirname "$targetprefix")"/temp_name.*
fi

exit 0