#!/bin/bash

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools
fi

# ----------------- Usage ------------------------------------------------

usage()
{
cat << EOF

Usage: subset -t <FILE> -a <FILE> [OPTIONS]

Subsets genome assembly file based on Blobtools database.

    -t  FILE            BlobtoolsDB file in table text format.
    -a  FILE            Genome assembly file to subset.

OPTIONS (Max. options at a time):
    -g  STRING          Group name to filter on ex. 'arthropoda'.
    -e  STRING          Exclude specific group.
    -m  RANGE           Range of GC content to filter on in as 'x-y'.
    -c  RANGE | MIN     Coverage to filter on either a range in format
                        'x-y' or a min. value in format 'x'.
    -l  RANGE | MIN     Sequence lenght to filter on either a range in
                        format 'x-y' or a min. value in format 'x'.
    
EOF
}

# ----------------- Flag Processing --------------------------------------

# if [ -z "$1" ]; then
#     usage
#     exit 1
# fi

param=0

while getopts 't:a:g:e:m:c:l:' OPTION; do
    case "$OPTION" in
		t)
	    	if [ -f "$OPTARG" ]; then
                targetdb="$(readlink -f "$OPTARG")"
                targetdir="$(dirname "$targetdb")"
                dbname=${targetdb%.*}
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
				exit 1 
            fi
	    	;;
        a)
            if [ -f "$OPTARG" ]; then
                targetasm="$(readlink -f "$OPTARG")"  
            else
                echo -e "\n$OPTARG does not seem to be a file\n"
				exit 1
            fi
            ;;
        g)
            group="$OPTARG"
            ((param++))
            ;;
        e)
            exgroup="$OPTARG"
            ((param++))
            ;;
        m)
            gcrange="$OPTARG"
            ((param++))
            ;;
        c)
            covrange="$OPTARG"
            ((param++))
            ;;
        l)
            lengthrange="$OPTARG"
            ((param++))
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done

if [ $((param)) -gt 1 ]; then
    echo -e "\nCan only filter on one parameter at a time (-g | -m | -c | -l)\n"
    exit 1
fi

# ----------------- Functions --------------------------------------------

asm_filter()
{
    # Name holder
    names=("" "no-" "gcrange_" "covrange_" "covmin_" "seqlenrange_" "seqlenmin_")

    awk '{OFS="\t"} {print $1}' "$dbname"_"${names[$1]}""$2".txt > "$targetdir"/temp_"${names[$1]}""$2"_idlist.txt

    count=
    while read -r line; do
        if [ -z $count ]; then
            rg -j 6 -F -A 1 "$line" "$targetasm" > "$targetdir"/"$speciesabbr"_"${names[$1]}""$2".fa
            count="1"
        else
            rg -j 6 -F -A 1 "$line" "$targetasm" >> "$targetdir"/"$speciesabbr"_"${names[$1]}""$2".fa
        fi
    done < "$targetdir"/temp_"${names[$1]}""$2"_idlist.txt

    rm -f "$targetdir"/temp_"${names[$1]}""$2"_idlist.txt
}

# ----------------- Filtering --------------------------------------------

speciesabbr="${dbname::6}"

# Filtering on group
if [ -n "$group" ]; then
    rg -j 6 '^# ' "$targetdb" > "$dbname"_"$group".txt
    rg -j 6 -F "${group^}" "$targetdb" >> "$dbname"_"$group".txt

    asm_filter 0 "$group" 
fi

# Excluding group
if [ -n "$exgroup" ]; then
    rg -j 6 '^# ' "$targetdb" > "$dbname"_no-"$exgroup".txt
    rg -j 6 -F -v "${exgroup^}" "$targetdb" >> "$dbname"_no-"$exgroup".txt

    asm_filter 1 "$exgroup"
fi

# Filtering on GC content
if [ -n "$gcrange" ]; then
    lower=${gcrange%-*}
    upper=${gcrange#*-}
    rg -j 6 '^# ' "$targetdb" > "$dbname"_gcrange_"$gcrange".txt
    awk -v lower="$lower" -v upper="$upper" 'BEGIN{OFS="\t" NR>12} ($3 >= lower && $3 <= upper) {print}' "$targetdb" >> "$dbname"_gcrange_"$gcrange".txt

    asm_filter 2 "$gcrange" 
fi

# Filtering on coverage
if [ -n "$covrange" ]; then
    # Coverage range
    if [[ "$covrange" == *-* ]]; then 
        lower=${covrange%-*}
        upper=${covrange#*-}
        rg -j 6 '^# ' "$targetdb" > "$dbname"_covrange_"$covrange".txt
        awk -v lower="$lower" -v upper="$upper" 'BEGIN{OFS="\t" NR>12} ($5 >= lower && $5 <= upper) {print}' "$targetdb" >> "$dbname"_covrange_"$covrange".txt

        asm_filter 3 "$covrange"

    # Coverage minimum
    else
        lower=$covrange
        rg -j 6 '^# ' "$targetdb" > "$dbname"_covmin_"$covrange".txt
        awk -v lower="$lower" 'BEGIN{OFS="\t" NR>12} ($5 >= lower) {print}' "$targetdb" >> "$dbname"_covmin_"$covrange".txt

        asm_filter 4 "$covrange"
    fi
fi

# Filtering on sequence length
if [ -n "$lengthrange" ]; then
    # Range of sequence lengths
    if [[ "$lengthrange" == *-* ]]; then 
        lower=${lengthrange%-*}
        upper=${lengthrange#*-}
        rg -j 6 '^# ' "$targetdb" > "$dbname"_seqlenrange_"$lengthrange".txt
        awk -v lower="$lower" -v upper="$upper" 'BEGIN{OFS="\t" NR>12} ($2 >= lower && $2 <= upper) {print}' "$targetdb" >> "$dbname"_seqlenrange_"$lengthrange".txt

        asm_filter 5 "$lengthrange" 

    # Minimum sequence length
    else
        lower=$lengthrange
        rg -j 6 '^# ' "$targetdb" > "$dbname"_covmin_"$lengthrange".txt
        awk -v lower="$lower" 'BEGIN{OFS="\t" NR>12} ($2 >= lower) {print}' "$targetdb" >> "$dbname"_seqlenmin_"$lengthrange".txt

        asm_filter 6 "$lengthrange"
    fi
fi

exit 0