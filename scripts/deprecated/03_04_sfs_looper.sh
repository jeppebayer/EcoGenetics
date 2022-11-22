#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

# ----------------- Queue Functions --------------------------------------

# Function to queue full sfs jobs
sfs_full()
{
    jid5=$(sbatch \
        --parsable \
        --time="$(timer "$9" 1800 120)" \ # ***Figure out time***
        --mem-per-cpu="$8"G \
        --cpus-per-task="$1" \
        --output="$stdoutput"/"$(basename "$5")"_05_"${13}"-%j.out \ # ***Figrue out to cycle names***
        "${10}"/modules/03_05_01_sfs_full.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7") # ***Add actual arguments***

    echo -e "\t'Site Frequency Spectrum - Full-${13}' job has been submitted for $(basename "$5") -- Job ID: $jid5" >> "$logfile" # ***Number to be passed as argument***
}

# Funtion to queue intergenic sfs jobs
sfs_intergenic()
{
    jid6=$(sbatch \
        --parsable \
        --time="$(timer "$9" 1800 120)" \ # ***Figure out time***
        --mem-per-cpu="$8"G \
        --cpus-per-task="$1" \
        --output="$stdoutput"/"$(basename "$5")"_06_"${13}"-%j.out \ # ***Figrue out to cycle names***
        "${10}"/modules/03_06_01_sfs_intergenic.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7") # ***Add actual arguments***

    echo -e "\t'Site Frequency Spectrum - Intergenic-${13}' job has been submitted for $(basename "$5") -- Job ID: $jid6" >> "$logfile" # ***Number to be passed as argument***
}

# Function to queue nonsynonymous sfs jobs
sfs_nonsyn()
{
    jid7=$(sbatch \
        --parsable \
        --time="$(timer "$9" 1800 120)" \ # ***Figure out time***
        --mem-per-cpu="$8"G \
        --cpus-per-task="$1" \
        --output="$stdoutput"/"$(basename "$5")"_07_"${13}"-%j.out \ # ***Figrue out to cycle names***
        "${10}"/modules/03_07_01_sfs_nonsyn.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7") # ***Add actual arguments***

    echo -e "\t'Site Frequency Spectrum - Non-synonymous-${13}' job has been submitted for $(basename "$5") -- Job ID: $jid7" >> "$logfile" # ***Number to be passed as argument***
}

# ----------------- Script Queue -----------------------------------------

# Counter to increment
n=1
# Looping through all parts to create full SFS
# "$cpus" "$RG" "$SD" "$WD" "$sample" "$dataprep" "$algo" "$memory" "$adjustment" "$script_path" "$age"
for part in "$WD"/temp/"$(basename "$sample")"_*.part; do
    sfs_full "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "$part" "$n"
    sfs_intergenic "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "$part" "$n"
    sfs_nonsyn "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "$part" "$n"
    ((n++))
done

exit 0