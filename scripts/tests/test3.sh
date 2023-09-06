#!/bin/bash

# Bed12 columns are 'chrom' |   'start' |   'end'   |   'name'  |   'score' |   'strand'    |   'thickStart'    |   'thickEnd'  |   'itemRgb'   |   'blockCount'    |   'blockSizes'    |   'blockStarts'
# From GTF              $1  |       $4  |   $5      |       $9  |       $6  |       $7      |       $4          |       $5      |       '0'     |       'N exons'   |   'diff start/end'|   'relative start pos'
awk -F '\t' '$3 == "transcript" && $1 != "#" {printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $4, $5, $9, $6, $7, $4, $5, "0")}' \
"${gtf_file[*]}" \
> "$tempcdsbed"

for ((i=1; i<=$(wc -l < "$tempcdsbed"); i++)); do
    column4=$(awk -F '\t' -v i=$i 'NR == i {print $4}' "$tempcdsbed")
    transcriptid=${column4#*transcript_id\ }
    transcriptid=${transcriptid%%;*}
    transcriptid=${transcriptid#*\"}
    transcriptid=${transcriptid%%\"*}

    cdschunks=("$(awk -F '\t' -v transcriptid="$transcriptid" \
    '$3 == "exon" && index($9, transcriptid) {
        ++n;
        lengtharray[n] = ($5 - $4);
        posarray[n] = $4
        } 
    END {
        j = lengtharray[1]
        sep = ","
        for (i=2; i<=n; i++) {
            j = j sep lengtharray[i]
            };
        k = (posarray[1] - posarray[1])
        for(i=2; i<=n; i++) {
            k = k sep (posarray[i] - posarray[1])
            };
        printf ("%s\t%s\t%s\n", n, j, k)
    }' \
    "${gtf_file[*]}")")

    awk -F '\t' -v transcriptid="$transcriptid" -v i=$i -v cdschunk1="${cdschunks[0]}" -v cdschunk2="${cdschunks[1]}" -v cdschunk3="${cdschunks[2]}" \
    'NR == i {
        printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, transcriptid, $5, $6, $7, $8, $9, cdschunk1, cdschunk2, cdschunk3)
        }' "$tempcdsbed" >> "$cdsbed12"
done

bedtools getfasta -fi "${reference_genome[*]}" -bed "$cdsbed" -nameOnly -split | fold -w 80 > "$(dirname "$cdsbed")"/cds.fa