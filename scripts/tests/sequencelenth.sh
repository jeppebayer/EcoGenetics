filepath=

for file in $filepath; do
    bioawk \
        -v OFS='\t' \
        -c fastx \
        '{ print length($seq) }' \
        <(zcat "$file") \
    | awk \
        '{Total=Total+$1}
        END{print "Total assembly length is: " Total}' \
        > "$(dirname "$(dirname "$file")")"/"$(basename "$file")".length
done