#!/bin/bash

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 05_repeatmm_00_init.sh -a <FILE> [OPTIONS]

    -a  FILE        Genome assembly.

OPTIONS:
    -d  PATH        Path to working directory. If none is provided the
                    directory of the referene genome will be used.
    -j  STRING      If 'list' shows jobgraph. Can be supplied with a single
                    integer or a range of integer to define which jobs in the
                    pipeline should be run.
    -h              Show this message.

EOF
}

# ----------------- Configuration ----------------------------------------

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command={print $2}')
# If run on the frontend:
else
    script=$(realpath "$0")
fi
script_path="$(dirname "$script")"

# Creates list for the execution of jobs
num_of_jobs=6
for (( i=1; i<="$num_of_jobs"; i++ )); do
    joblist["$i"]=true
done

jobgraph()
{
cat << EOF

Sequential dependence:
All previous steps must have been run but sequence can be stopped after
any particular step.

        1   Create RepeatModeler DB
                        | |
        2   Run RepeatModeler on Assembly
                        | |
        3   Run RepeatMasker on Assembly using RepBase DB
                        | |
        4   Run RepeatMasker on Assembly using RepeatModeler data
                        | |
        5   Combine result of 3 and 4 and Process them to create new files
                        | |
        6   Mask Assembly based on 5

EOF
}

# ------------------ Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 1
fi

while getopts 'a:d:j:h' OPTION; do
    case "$OPTION" in
		a)
	    	if [ -f "$OPTARG" ]; then
                reference_genome="$(readlink -f "$OPTARG")"
                referencedir="$(dirname "$reference_genome")"
                species_name="$(basename "$referencedir")"
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
				exit 1 
            fi
	    	;;
        d)
            WD="$(readlink -f "$OPTARG")"  
            ;;
        j)
            if [ "$OPTARG" == "list" ]; then
                jobgraph
                exit 0
            else
                for (( i=1; i<="$num_of_jobs"; i++ )); do
                    joblist["$i"]=false
                done
                if [ ${#OPTARG} -gt 2 ]; then
                    min=${OPTARG%-*}
                    max=${OPTARG#*-}
                    for (( i="$min"; i<="$max"; i++ )); do
                        joblist["$i"]=true
                    done
                else
                    joblist["$OPTARG"]=true
                fi
            fi
            ;;
        h)
            usage
            exit 0
            ;;
        ?)
            usage
            exit 2
            ;;
    esac
done

# Tests that reference genome has been supplied
if [[ -z "$reference_genome" ]]; then
    usage
    exit 3
fi

# If a working directory is not supplied the directory of the reference genome will be used
if [ -z "$WD" ]; then
    WD="$referencedir"
else
    # Removes "/" from the end of path for working directory
    if [ "${WD: -1}" == "/" ]; then
        length=${#WD}
        WD=${WD::length - 1}
    fi
    # Creates working directory if it doesn't exist
    [ -d "$WD" ] || mkdir -m 775 "$WD"
    # Creates new working folder
    WD="$WD"/05_repeatmm
    [ -d "$WD" ] || mkdir -m 775 "$WD"
    # Creates new species folder
    WD="$WD"/"$species_name"
    [ -d "$WD" ] || mkdir -m 775 "$WD"
fi

# ----------------- Functions --------------------------------------------

# Function to flexibly control the dependencies of each job
dependency(){
    depend_on=
    for id in $1; do
         if [ -n "${jobidlist[$id]}" ]; then
            if [ -z "$depend_on" ]; then
                depend_on="${jobidlist[$id]}"
            else
                depend_on="$depend_on:${jobidlist[$id]}"
            fi
        fi
    done
    [ -z "$depend_on" ] && depend_on="-1"
    echo "$depend_on"
}

job1(){
    jobidlist[1]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$rmdb" \
                --time=30 \
                --mem=16G \
                --cpus-per-task=2 \
                --dependency=afterany:"$(dependency "")" \
                --output="$out"/01_database-%j.out \
                "$script_path"/modules/05_repeatmm_01_database.sh "$reference_genome")
}

job2(){
    # Will take a long time and a lot of resources
    jobidlist[2]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$repmod" \
                --time=4320 \
                --mem=192G \
                --cpus-per-task=32 \
                --dependency=afterany:"$(dependency "1")" \
                --output="$out"/02_model-%j.out \
                "$script_path"/modules/05_repeatmm_02_model.sh "$reference_genome")
}

job3(){
    jobidlist[3]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$repbaserun" \
                --time=300 \
                --mem=192G \
                --cpus-per-task=32 \
                --dependency=afterany:"$(dependency "2")" \
                --output="$out"/03_mask_repbase-%j.out \
                "$script_path"/modules/05_repeatmm_03_mask_repbase.sh "$reference_genome" "$repbaserun" "$repmodrun" "$repmod")
}

job4(){
    jobidlist[4]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$repmodrun" \
                --time=300 \
                --mem=192G \
                --cpus-per-task=32 \
                --dependency=afterany:"$(dependency "3")" \
                --output="$out"/04_mask_repmod-%j.out \
                "$script_path"/modules/05_repeatmm_04_mask_repmod.sh "$reference_genome" "$repbaserun" "$repmodrun" "$repmod")
}

job5(){
    jobidlist[5]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$repmask" \
                --time=60 \
                --mem=32G \
                --cpus-per-task=4 \
                --dependency=afterany:"$(dependency "4")" \
                --output="$out"/05_combine_process-%j.out \
                "$script_path"/modules/05_repeatmm_05_combine_process.sh "$reference_genome" "$repbaserun" "$repmodrun")
}

job6(){
    jobidlist[6]=$(sbatch \
                --parsable \
                --account=EcoGenetics \
                --chdir="$repmask" \
                --time=60 \
                --mem=32G \
                --cpus-per-task=4 \
                --dependency=afterany:"$(dependency "5")" \
                --output="$out"/06_maskfasta-%j.out \
                "$script_path"/modules/05_repeatmm_06_maskfasta.sh "$reference_genome" "$repmask" "$WD")
}

# ----------------- Repeat Masking ---------------------------------------

# Creates directory for standard out files
out=$WD/out
[ -d "$out" ] || mkdir -m 775 "$out"

# Creates directory for RepeatModeler
repmod="$WD"/RepeatModeler
[ -d "$repmod" ] || mkdir -m 775 "$repmod"

# Creates directory for RepeatModeler database
rmdb="$repmod"/RM_DB_"$species_name"
[ -d "$rmdb" ] || mkdir -m 775 "$rmdb"

# Creates directory for RepeatMasker
repmask="$WD"/RepeatMasker
[ -d "$repmask" ] || mkdir -m 775 "$repmask"

# Creates directory for RepeatMasker run on RepBase Data
repbaserun="$repmask"/RepBase_run
[ -d "$repbaserun" ] || mkdir -m 775 "$repbaserun"

# Creates directory for RepeatMasker run on RepModeler Data
repmodrun="$repmask"/RepMod_run
[ -d "$repmodrun" ] || mkdir -m 775 "$repmodrun"

logmessage=""

# Executes job if indicated in joblist
[ "${joblist[1]}" = true ] && job1 && logmessage+="1\tCreate RepeatModeler DB\tJobID:\t${jobidlist[1]}\n"

# Executes job if indicated in joblist
[ "${joblist[2]}" = true ] && job2 && logmessage+="2\tRun RepeatModeler on Assembly\tJobID:\t${jobidlist[2]}\n"

# Executes job if indicated in joblist
[ "${joblist[3]}" = true ] && job3 && logmessage+="3\tRun RepeatMasker on Assembly using RepBase DB\tJobID:\t${jobidlist[3]}\n"

# Executes job if indicated in joblist
[ "${joblist[4]}" = true ] && job4 && logmessage+="4\tRun RepeatMasker on Assembly using RepeatModeler data\tJobID:\t${jobidlist[4]}\n"

# Executes job if indicated in joblist
[ "${joblist[5]}" = true ] && job5 && logmessage+="5\tCombine result of 3 and 4 and Process them to create new files\tJobID:\t${jobidlist[5]}\n"

# Executes job if indicated in joblist
[ "${joblist[6]}" = true ] && job6 && logmessage+="6\tMask Assembly based on 5\tJobID:\t${jobidlist[6]}\n"

log="$WD"/repeatmm_"$species_name".log
echo -n "" > "$log"

echo -e "Start date:\t$(date +"%Y-%m-%d")\nTarget assembly:\t$reference_genome\nJobs sent to queue are as follows:\n\n$logmessage" >> "$log"

echo -e "\nLogfile can be found here: $log\n"

exit 0