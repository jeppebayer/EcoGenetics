#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:40:00
#SBATCH --output=Index_Reference_Genome-%j.out

if [ "$USER" == "jepe" ]; then

    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary

fi

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: 01_initialize_indexing_reference.sh [OPTIONS]

This script is used for initializing indexation of reference genomes used at
the Center for EcoGenetics.
Can be run on the frontend, as its resource-demands are fairly, or run in 
conjunction with 'sbatch'or 'srun'.

OPTIONS:
    -r  FILE            Reference genome file
    -f                  Force run. Index reference genome even if it has been
                        indexed previously
    -a                  Index all available reference genomes (Cannot be run with -r)
    -h                  Show this message

Tested and working using:
'EcoGenetics/people/Jeppe_Bayer/environment_primary_from_history.yml'

EOF
}

# ----------------- Configuration ----------------------------------------

# Option to force run even if target sample directory contains index files
force_overwrite="N"

# Option to run on all available reference genomes
run_all="N"

# Gets path to script
# If run through sbatch or srun:
if [ -n "$SLURM_JOB_ID" ]; then
    script_path=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
# If run on the frontend:
else
    script_path=$(realpath "$0")
fi

# ----------------- Script Flag Processing -------------------------------

while getopts 'r:fah' OPTION; do
    case "$OPTION" in
        r)
            if [ -f "$OPTARG" ]; then
                RG="$(readlink -f "$OPTARG")"
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
                exit 1
            fi
            ;;
        f)
            force_overwrite="Y"
            ;;
        a)
            run_all="Y"
            ;;
        h)
            usage
            exit 1
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done

# # Removes "/" from the end of path for reference genome directory
# if [ "${RGD: -1}" == "/" ]; then
#     length=${#RGD}
#     RG=${RGD::length - 1}
# fi

# ----------------- Functions --------------------------------------------
queue()
{
    # Index Reference
	jid1=$(sbatch \
		--parsable \
		--time=120 \
		--mem-per-cpu=16G \
		--cpus-per-task=2 \
		--output=/dev/null \
		--error=/dev/null \
		"$script_path"/modules/01_01_indexing_reference.sh "$RG")
}

# ----------------- Script Queue -----------------------------------------

# Holds path to script directory
script_path=$(dirname "$script_path")

# Makes sure all associated modules are executable
for script in "$script_path"/modules/*; do
    [ ! -x "$script" ] && chmod -R u+rwx "$script_path"/modules && break
done

# Check whether or not to run on all available reference genomes
if [ "$run_all" == "N" ]; then
    
    if [[ "$RG" == *.fna ]] || [[ "$RG" == *.fasta ]] || [[ "$RG" == *.fa ]]; then

        # Exit with error message if reference genome is indexed
        for index in "$(dirname "$RG")"/*.ann; do
            
            if [ -e "$index" ] && [ "$force_overwrite" == "N" ]; then

                echo -e "\nDesignated reference genome, $(basename "$RGD"), is already indexed\n"
                exit 1

            else
                
                queue
                break

            fi

        done

    else

        echo -e "\nDesignated reference genome, $(basename "$RG"), is not in .fna, .fa or .fasta format\n"
        exit 1
        
    fi

else

    for RGD in /faststorage/project/EcoGenetics/BACKUP/reference_genomes/*; do

        # Checks that folder contains .fna file
        for fna in "$RGD"/*.fna; do

            if [ -e "$fna" ]; then

                # Skip if reference genome is indexed
                for index in "$RGD"/*.ann; do

                    if [ -e "$index" ] && [ "$force_overwrite" == "N" ]; then

                        echo -e "\nReference genome, $(basename "$RGD"), is already indexed and is skipped"
                        break

                    else

                        for RG in "$RGD"/*.fna; do

                            queue
                            break

                        done

                    break

                    fi

                done

            else

                echo -e "\n$RGD does not contain a .fna file and is skipped"
                break

            fi  

        done

    done

fi

exit 0