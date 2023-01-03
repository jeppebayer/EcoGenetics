#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 2
#SBATCH --time 00:30:00
#SBATCH --output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/04_assembly/gfa_to_fasta-%j.out

if [[ ! $1 == *.gfa ]]; then
    echo "$1 is not in .gfa format"
    exit 1
fi

gfa_in="$(readlink -f "$1")"
fasta_out="${gfa_in%.*}"

awk '/^S/{print ">"$2"\n"$3}' "$gfa_in" | fold > "$fasta_out".fasta

exit 0