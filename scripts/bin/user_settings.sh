#!/bin/bash

test_path="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/bin"

if [ -f "$HOME"/.bash_profile ]; then
    # echo -e "export PATH=$PATH:/faststorage/project/EcoGenetics/scripts" >> "$HOME"/.bash_profile
    echo ".bash_profile found... Adding path..."
    echo -e "export PATH=\$PATH:$test_path" >> "$HOME"/.bashrc
else
    # echo -e "# .bash_profile\n\n# Get the aliases and functions\nif [ -f ~/.bashrc ]; then\n\t. ~/.bashrc\nfi\n\n# User specific environment and startup programs\n\nexport PATH=$PATH:/faststorage/project/EcoGenetics/scripts" > "$HOME"/.bash_profile
    echo ".bash_profile doesn't exist... Creating..."
    echo -e "# .bash_profile\n\n# Get the aliases and functions\nif [ -f ~/.bashrc ]; then\n\t. ~/.bashrc\nfi\n\n# User specific environment and startup programs\n\nexport PATH=\$PATH:$test_path" > "$HOME"/.bash_profile
fi

source "$HOME"/.bashrc