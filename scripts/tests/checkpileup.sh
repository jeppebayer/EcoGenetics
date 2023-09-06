#!/bin/bash

pileup_file=/faststorage/project/EcoGenetics/people/EliseL/Functional_neutral_GD/Populations/Outputs/Orchesella_cincta/NYS-F/orc_cinc_NYSF.pileup
pileup_name="$(basename "$pileup_file")"
pileup_name=${pileup_name%.*}

awk \
    -F "\t" \
    'BEGIN {
        OFS="\t" ; 
        print "Contig", "Position", "N_reads", "N_qualities", "Case"
    }
    {if (($4 == "0") && ($6 == "*"))
        {answer = "Del"}
    else 
        {if ($4 == length($6))
            {answer = "Match"}
        else
            {answer = "Mismatch"}
        }
    if (answer == "Mismatch")
        {print ($1, $2, $4, length($6), answer)
         count += 1}
    }
    END {print "Total", "number", "of", "mismatches", count}' \
    "$pileup_file" \
    > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/"$pileup_name".check.txt

exit 0