#!/bin/bash

reference_genome="$1"
repbaserun="$2"
repmodrun="$3"
repmod="$4"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate repeatmasking
fi

RepeatMasker \
-e rmblast \
-pa 8 \
-dir "$repbaserun" \
-xsmall \
-lib /faststorage/project/EcoGenetics/BACKUP/database/giri_repbase/RepBaseCustom15.02.fasta/Arthropoda.rep \
"$reference_genome"

# Make GFF 3 of repeat regions from RepeatMasker output file
repeatmasker_out="$repbaserun"/"$(basename "$reference_genome")".out
filename="$(basename "$repeatmasker_out")"
filename=${filename%.*.*}

awk \
    -F " " \
    'BEGIN {OFS="\t"; print "##gff-version 3"}
    NR > 3
        {if ($9 == "C")
            {strand = "-"}
        else
            {strand = "+"}
        if ($12 ~ /\(/)
            {start = $14}
        else {start = $12}
        print $5, "RepeatMasker", "repeat_region", $6, $7, ".", strand, ".", "ID="$15";Name="$10";Class="$11";Family="$11";Target="$10" "start" "$13}' \
    "$repeatmasker_out" \
    > "$(dirname "$repeatmasker_out")"/"$filename".repeats_giri.gff

exit 0

# -lib 
# 	The majority of species are of course not yet covered in the repeat
# 	databases and many are far from complete, but you may have your own
# 	collection. At other times you may want to mask or study only a
# 	particular type of repeat.

# 	For these types of siutations, you can use the -lib option to
# 	specify a custom library of sequences to be masked in the query. The
# 	library file needs to contain sequences in FASTA format. Unless a full
# 	path is given on the command line the file is assumed to be in the
# 	same directory as the sequence file.  

# 	The recommended format for IDs in a custom library is:

# 	>repeatname#class/subclass
# 	or simply
# 	>repeatname#class

# 	In this format, the data will be processed (overlapping repeats are
# 	merged etc), alternative output (.ace or .gff) can be created and an
# 	overview .tbl file will be created. Classes that will be displayed in
# 	the .tbl file are 'SINE', 'LINE', 'LTR', 'DNA', 'Satellite', anything
# 	with 'RNA' in it, 'Simple_repeat', and 'Other' or 'Unknown' (the
# 	latter defaults when class is missing). Subclasses are plentiful. They
# 	are not all tabulated in the .tbl file or necessarily spelled
# 	identically as in the repeat files, so check the RepeatMasker.embl
# 	file for names that can be parsed into the .tbl file.

# 	You can combine the repeats available in the RepeatMasker library 
# 	with a custom set of consensus sequences.  To accomplish this 
# 	use the famdb.py tool:

# 	`./famdb.py -i Libraries/RepeatMaskerLib.h5 families --format fasta_name --ancestors --descendants 'species name' --include-class-in-name`

# 	The resulting sequences can be concatenated to your own set of sequences in a
# 	new library file.

# -cutoff
# 	When using a local library you may want to change the minimum score
# 	for reporting a match. The default is 225, lowering it below 200 will
# 	usually start to give you significant numbers of false matches,
# 	raising it to 250 will guarantee that all matches are real. Note that
# 	low complexity regions in otherwise complex repeat sequences in your
# 	library are most likely to give false matches.

# -small
# 	returns complete .masked sequence in lower case

# -xsmall
# 	When the option -xsmall is used a sequence is returned in the .masked
# 	file in which repeat regions are in lower case and non-repetitive
# 	regions are in capitals.

# -gff
# 	The script creates a .gff file with the annotation in 'General Feature
# 	Finding' format. See http://www.sanger.ac.uk/Software/GFF for
# 	details. The current output follows a Sanger convention:

# 	<seqname> RepeatMasker Similarity <start in query> <end in query>
# 	<percent divergence> <orientation> . Target "Motif:<repeat-name>"
# 	<start in consensus> <end in consensus>

# 2.3  Repeat databases

# The RepeatMasker program are distributed with a copy of the
# Dfam database ( www.dfam.org ). Dfam is a small but growing "open"
# databases of Transposable Element seed alignments, profile Hidden
# Markov Models and consensus sequences.

# RepeatMasker is also compatible with the RepBase database managed by
# the Genetic Information Research Institute and requires a license to
# use. Up until 2019 we maintained the "Repbase RepeatMasker Edition"
# libraries as co-editor of RepBase Update.  For newer versions of
# RepBase users will need to use the sequences in FASTA format with
# RepeatMasker's "-lib" option.

# 2.4 Sensitivity and speed  

# Note that for many non-mammalian species the slower settings do not
# dramatically increase the percentage recognized as interspersed
# repeats. Most of the repeats in the databases for these species are
# relatively young and thus are easily detected. This particular 1Mbp
# Arabidopsis sequence is an extreme example, where at slow settings in
# almost two hours only 1800 bp more is masked than at rush settings in
# 51 seconds (the Arabidopsis database is large).

# The user time for larger sequences or sequence batches (50 kb and up)
# is linearly related to the length of the query due to the
# fragmentation of the query sequence.