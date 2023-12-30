#!/bin/bash

# ------------------ Usage ------------------------------------------------

usage(){
cat << EOF

Usage: $(basename "$0") [-n <PROJECT_NAME>]

Setup directory structure for new project.

OPTIONS:
    -n  STRING          Desired name of project directory.
                        Project name can contain path.
                        Default: 'project'.
    -h                  Show this message.

EOF
}

# ------------------ Configuration ----------------------------------------

project_name=project
directory_path="$PWD"

# ------------------ Functions --------------------------------------------

make_directories(){
    echo -en "\nDo you want to create the project directory '$(basename "$project_name")' in the directory\n'$directory_path'?"
    read -p " (Y/n): " answer
    if [[ "$answer" == @("yes"|"Yes"|"y"|"Y"|"") ]]; then
        mkdir -m 764 "$directory_path"/"$project_name"
        mkdir -m 774 "$directory_path"/"$project_name"/scripts
        mkdir -m 764 "$directory_path"/"$project_name"/steps
        mkdir -m 764 "$directory_path"/"$project_name"/results
        echo "Project directory has been created."
        exit 0
    else
        echo "Cancelling..."
        exit 0
    fi
}

# ------------------ Flag Processing --------------------------------------

if [ -z "$1" ]; then
    make_directories
fi

if [ "$#" -gt 2 ]; then
    echo "Too many arguments"
    exit 3
fi

while getopts 'n:h' OPTION; do
    case "$OPTION" in
        
        n)
            project_name="$OPTARG"
            project_name="${project_name/" "/"_"}"
            if [[ "$project_name" == */* ]]; then
                directory_path=$(readlink -f "$(dirname "$project_name")")
                if [ -d "$directory_path" ]; then
                    make_directories
                else
                    echo -e "\n$directory_path does not seem to exist\n"
                    exit 1
                fi
            else
                make_directories
            fi
            ;;
        h)
            usage
            exit 0
            ;;
        ?)
            echo -e "\nUnrecognized option"
            usage
            exit 2
            ;;
    esac
done

echo -e "\n$(basename "$0") does not accept positional arguments."
usage
exit 4