#!/bin/bash

bedfile="$(readlink -f "$1")"
filename="$(basename "$bedfile")"
filename=${filename%.*}

awk \
    -F "\t" \
    'BEGIN {OFS="\t"; print "##gtf-version 2.2"}
    {print $1, "Custom_script", "inter", ($2 + 1), $3, ".", ".", ".", "gene_id \"\"; transcript_id \"\";"}' \
    "$bedfile" \
    > "$(dirname "$bedfile")"/"$filename".gtf

exit 0