#!/bin/bash

gff_file=/faststorage/project/EcoGenetics/BACKUP/reference_genomes/Aglais_urticae/annotation/Nymphalis_urticae-GCA_905147175.2-2022_03-genes.gff3
output=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/gff_name_list.txt
temp_gff=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/tests/gff_temp.gff
new_gff=/faststorage/project/EcoGenetics/BACKUP/reference_genomes/Aglais_urticae/annotation/Nymphalis_urticae-GCA_905147175.2-2022_03-genes.rename.gff3

# grep \
#     'Alias=' \
#     "$gff_file" \
# | awk \
#     -F '[\t=]' \
#     '{print $1 "\t" $11}' \
#     > $output

awk \
    -F "\t" \
    'BEGIN{OFS="\t"}
    {if ($1 == 1)
        {gsub("1", "LR989983.1", $1)}
    if ($1 == 10)
        {gsub("10", "LR989992.1", $1)}
    if ($1 == 11)
        {gsub("11", "LR989993.1", $1)}
    if ($1 == 12)
        {gsub("12", "LR989994.1", $1)}
    if ($1 == 13)
        {gsub("13", "LR989995.1", $1)}
    if ($1 == 14)
        {gsub("14", "LR989997.1", $1)}
    if ($1 == 15)
        {gsub("15", "LR989998.1", $1)}
    if ($1 == 16)
        {gsub("16", "LR989999.1", $1)}
    if ($1 == 17)
        {gsub("17", "LR990000.1", $1)}
    if ($1 == 18)
        {gsub("18", "LR990001.1", $1)}
    if ($1 == 19)
        {gsub("19", "LR990002.1", $1)}
    if ($1 == 2)
        {gsub("2", "LR989984.1", $1)}
    if ($1 == 20)
        {gsub("20", "LR990003.1", $1)}
    if ($1 == 21)
        {gsub("21", "LR990004.1", $1)}
    if ($1 == 22)
        {gsub("22", "LR990005.1", $1)}
    if ($1 == 23)
        {gsub("23", "LR990006.1", $1)}
    if ($1 == 24)
        {gsub("24", "LR990007.1", $1)}
    if ($1 == 25)
        {gsub("25", "LR990008.1", $1)}
    if ($1 == 26)
        {gsub("26", "LR990009.1", $1)}
    if ($1 == 27)
        {gsub("27", "LR990010.1", $1)}
    if ($1 == 28)
        {gsub("28", "LR990011.1", $1)}
    if ($1 == 29)
        {gsub("29", "LR990012.1", $1)}
    if ($1 == 3)
        {gsub("3", "LR989985.1", $1)}
    if ($1 == 30)
        {gsub("30", "LR990013.1", $1)}
    if ($1 == 4)
        {gsub("4", "LR989986.1", $1)}
    if ($1 == 5)
        {gsub("5", "LR989987.1", $1)}
    if ($1 == 6)
        {gsub("6", "LR989988.1", $1)}
    if ($1 == 7)
        {gsub("7", "LR989989.1", $1)}
    if ($1 == 8)
        {gsub("8", "LR989990.1", $1)}
    if ($1 == 9)
        {gsub("9", "LR989991.1", $1)}
    if ($1 == "CAJHUP020000001.1")
        {gsub("CAJHUP020000001.1", "SUPER_26_unloc_1", $1)}
    if ($1 == "CAJHUP020000002.1")
        {gsub("CAJHUP020000002.1", "SUPER_26_unloc_2", $1)}
    if ($1 == "W")
        {gsub("W", "LR989996.1", $1)}
    if ($1 == "Z")
        {gsub("Z", "LR989982.1", $1)}
    print $0}' \
    $gff_file \
    > $new_gff
