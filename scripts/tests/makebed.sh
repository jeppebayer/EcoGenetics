#!/bin/bash

# ----------------- Configuration ----------------------------------------

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate vcf
fi

# ----------------- Usage ------------------------------------------------

usage(){
cat << EOF

Usage: makebed.sh <OPTION>

OPTIONS:
    intergenic_bed2gtf                  Create GTF 2.2 file of intergenic
                                        regions in BED file
    gene_intergene_gff2bed              Create two BED files, one of gene
                                        regions and one of intergenic
                                        regions from GFF file
    repeatout2gff                       Make GFF 3 of repeat regions from
                                        RepeatMasker output file
    gff2bed                             Make BED file of GFF 3 file

EOF
}

# ------------------ Flag Processing --------------------------------------

if [ -z "$1" ]; then
    usage
    exit 0
fi

case "$1" in
    intergenic_bed2gtf)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: makebed.sh intergenic_bed2gtf <BED_FILE>\n\nCreate GTF 2.2 file of intergenic regions in BED file\n"
        else
            # Create GTF 2.2 file of intergenic region in BED file
            bed_file="$(readlink -f "$2")"
            filename="$(basename "$bed_file")"
            filename=${filename%.*}

            awk \
                -F "\t" \
                'BEGIN {OFS="\t"; print "##gtf-version 2.2"}
                {print $1, "Custom_script", "inter", ($2 + 1), $3, ".", ".", ".", "gene_id \"\"; transcript_id \"\";"}' \
                "$bed_file" \
                > "$(dirname "$bed_file")"/"$filename".gtf
        fi
        ;;
    gene_intergene_gff2bed)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: makebed.sh gene_intergene_gff2bed <GFF_FILE> <REFERENCE_GENOME>\n\nCreate two BED files, one of gene regions and one of intergenic regions from GFF file\n"
        else
            # Create two BED files, one of gene regions and one of intergenic regions from GFF file.
            gff_file="$(readlink -f "$2")"
            genome_sequence="$(readlink -f "$3")"
            filename="$(basename "$gff_file")"
            filename=${filename%.*}

            grep \
                'gene' \
                "$gff_file" \
            | awk \
                -F "\t" \
                'BEGIN {OFS = "\t"}
                {if ($0 ~/^#/)
                    {next;}
                print $1, ($4 - 1), $5}' \
            | sort \
                -k1,1 \
                -k2,2n \
                > "$(dirname "$gff_file")"/"$filename".gene.bed

            bedtools complement \
                -i "$(dirname "$gff_file")"/"$filename".gene.bed \
                -g <(bioawk \
                        -c fastx \
                        'BEGIN {OFS="\t"}
                        {print $name, length($seq)}' \
                        "$genome_sequence" \
                    | sort \
                        -k1,1 \
                        -k2,2n) \
                > "$(dirname "$gff_file")"/"$filename".intergenic.bed
        fi
        ;;
    repeatout2gff)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: makebed.sh repeatout2bed <RepeatMasker.out>\n\nMake GFF 3 of repeat regions from RepeatMasker output file\n"
        else
            # Make GFF 3 of repeat regions from RepeatMasker output file
            repeatmasker_out="$(readlink -f "$2")"
            filename="$(basename "$repeatmasker_out")"
            filename=${filename%.*.*.*}
            
            awk \
                -F " " \
                'BEGIN {OFS="\t"; print "##gff-version 3"}
                NR > 3 {
                    if ($9 == "C")
                        {strand = "-"}
                    else
                        {strand = "+"}
                    if ($12 ~ /\(/)
                        {start = $14}
                    else {start = $12}
                    print $5, "RepeatMasker", "repeat_region", $6, $7, ".", strand, ".", "ID="$15";Name="$10";Class="$11";Family="$11";Target="$10" "start" "$13
                }' \
                "$repeatmasker_out" \
                > "$(dirname "$repeatmasker_out")"/"$filename".repeats.gff
        fi
        ;;
    gff2bed)
        if [ "$#" == 1 ]; then
            echo -e "\nUsage: makebed.sh gff2bed <GFF_FILE>\n\nMake BED file of GFF 3 file\n"
        else
            # Make BED file of GFF 3 file
            gff_file="$(readlink -f "$2")"
            filename="$(basename "$gff_file")"
            filename=${filename%.*}

            awk \
                -F "\t" \
                'BEGIN {OFS = "\t"}
                {if ($0 ~ /^[^#]/)
                    {print $1, ($4 - 1), $5}
                }' \
                "$gff_file" \
                > "$(dirname "$gff_file")"/"$filename".repeats.bed
        fi
        ;;
    h)
        usage
        exit 0
        ;;
    ?)
        usage
        exit 2
        ;;
esac

exit 0

# GFF 3
# |   1   |   2    |  3   |   4   |  5  |   6   |   7    |   8   |     9      |
# | seqid | source | type | start | end | score | strand | phase | attributes |
# 1 : The ID of the landmark used to establish the coordinate system for the current feature.
# 2 : The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature. Typically this is the name of a piece of software, such as "Genescan" or a database name, such as "Genbank."
# 3 : The SOFA (Sequence Ontology Feature Annotation) feature type most equivalent to the feature found in the source annotation.
# 4 : The start coordinates of the feature are given in positive 1-based integer coordinates, relative to the landmark given in column one. Start is always less than or equal to end. For features that cross the origin of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes), the requirement for start to be less than or equal to end is satisfied by making end = the position of the end + the length of the landmark feature. For zero-length features, such as insertion sites, start equals end and the implied site is to the right of the indicated base in the direction of the landmark.
# 5 : The end coordinates of the feature are given in positive 1-based integer coordinates, relative to the landmark given in column one.
# 6 : The score of the feature, a floating point number. As in earlier versions of the format, the semantics of the score are ill-defined. It is strongly recommended that E-values be used for sequence similarity features, and that P-values be used for ab initio gene prediction features.
# 7 : The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features that are not stranded. In addition, ? can be used for features whose strandedness is relevant, but unknown.
# 8 : For features of type "CDS", the phase indicates where the next codon begins relative to the 5' end (where the 5' end of the CDS is relative to the strand of the CDS feature) of the current CDS feature. For clarification the 5' end for CDS features on the plus strand is the feature's start and and the 5' end for CDS features on the minus strand is the feature's end. The phase is one of the integers 0, 1, or 2, indicating the number of bases forward from the start of the current CDS feature the next codon begins. A phase of "0" indicates that a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward), a phase of "1" indicates that the codon begins at the second nucleotide of this CDS feature and a phase of "2" indicates that the codon begins at the third nucleotide of this region. Note that ‘Phase’ in the context of a GFF3 CDS feature should not be confused with the similar concept of frame that is also a common concept in bioinformatics. Frame is generally calculated as a value for a given base relative to the start of the complete open reading frame (ORF) or the codon (e.g. modulo 3) while CDS phase describes the start of the next codon relative to a given CDS feature. The phase is REQUIRED for all CDS features.
# 9 : A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons. URL escaping rules are used for tags or values containing the following characters: ",=;". Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape. Attribute values do not need to be and should not be quoted. The quotes should be included as part of the value by parsers and not stripped.

# These tags have predefined meanings:
# ID: Indicates the ID of the feature. The ID attribute is required for features that have children (e.g. gene and mRNAs), or for those that span multiple lines, but are optional for other features. IDs for each feature must be unique within the scope of the GFF file. In the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines. All lines that share an ID must collectively represent a single feature.
# Name: Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file.
# Alias: A secondary name for the feature. It is suggested that this tag be used whenever a secondary identifier for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias be unique within the file.
# Parent: Indicates the parent of the feature. A parent ID can be used to group exons into transcripts, transcripts into genes, an so forth. A feature may have multiple parents. Parent can only be used to indicate a partof relationship.
# Target: Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. The format of the value is "target_id start end [strand]", where strand is optional and may be "+" or "-". If the target_id contains spaces, they must be escaped as hex escape %20.
# Gap: The alignment of the feature to the target if the two are not collinear (e.g. contain gaps). The alignment format is inspired from the CIGAR format described in the Exonerate documentation.
# Derives_from: Used to disambiguate the relationship between one feature and another when the relationship is a temporal one rather than a purely structural "part of" one. This is needed for polycistronic genes. See "PATHOLOGICAL CASES" for further discussion.
# Note: A free text note.
# Dbxref: A database cross reference. See the section "Ontology Associations and Db Cross References" for details on the format.
# Ontology_term: A cross reference to an ontology term. See the section "Ontology Associations and Db Cross References" for details.
# Is_circular: A flag to indicate whether a feature is circular. See extended discussion below. 

# GTF 2.2
# |    1    |   2    |    3    |   4   |  5  |   6   |   7    |   8   |      9       |     10     |
# | seqname | source | feature | start | end | score | strand | frame | *attributes* | *comments* |
# 1 : The name of the sequence. Commonly, this is the chromosome ID or contig ID
# 2 : The source column should be a unique label indicating where the annotations came from
# 3 : The following feature types are required: "CDS", "start_codon", "stop_codon". The features "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are optional. All other features will be ignored. The types must have the correct capitalization shown here
# 4 : Integer start coordinates of the feature relative to the beginning of the sequence named in <seqname>.  <start> must be less than or equal to <end>. Sequence numbering starts at 1. Values of <start> and <end> that extend outside the reference sequence are technically acceptable, but they are discouraged.
# 5 : Integer end coordinates of the feature relative to the beginning of the sequence named in <seqname>
# 6 : The score field indicates a degree of confidence in the feature's existence and coordinates. The value of this field has no global scale but may have relative significance when the <source> field indicates the prediction program used to create this annotation. It may be a floating point number or integer, and not necessary and may be replaced with a dot. 
# 7 : Which stand feature is located on: + or -
# 8 : 0 indicates that the feature begins with a whole codon at the 5' most base. 1 means that there is one extra base (the third base of a codon) before the first whole codon and 2 means that there are two extra bases (the second and third bases of the codon) before the first codon. Note that for reverse strand features, the 5' most base is the <end> coordinate.
# 9 : All nine features have the same two mandatory attributes at the end of the record: gene_id value; A globally unique identifier for the genomic locus of the transcript. If empty, no gene is associated with this feature. transcript_id value; A globally unique identifier for the predicted transcript. If empty, no transcript is associated with this feature. These attributes are designed for handling multiple transcripts from the same genomic region. Any other attributes or comments must appear after these two and will be ignored. Attributes must end in a semicolon which must then be separated from the start of any subsequent attribute by exactly one space character (NOT a tab character). Textual attributes should be surrounded by doublequotes. These attributes are required even for non-mRNA transcribed regions such as "inter" and "inter_CNS" features.
# 10 : Comments begin with a hash ('#') and continue to the end of the line. Nothing beyond a hash will be parsed. These may occur anywhere in the file, including at the end of a feature line.

# RepeatMasker Output
# 1 : Smith-Waterman score of the match
# 2 : % substitutions in matching region compared to the concensus
# 3 : % of bases opposite a gap in the query sequence (deleted bp)
# 4 : % of bases opposite a gap in the repeat consensus (inserted bp)
# 5 : Name of query sequence
# 6 : Starting position of match in query sequence (1-based)
# 7 : Ending position of match in query sequence
# 8 : No. of bases in query sequence past the ending position of match
# 9 : Match is with the Complement of the consensus sequence in the database
# 10 : Name of the matching interspersed repeat
# 11 : The class of the repeat
# 12 : No. of bases in (complement of) the repeat consensus sequence prior to beginning of the match (so 0 means that the match extended all the way to the end of the repeat consensus sequence)
# 13 : Starting position of match in database sequence (using top-strand numbering)
# 14 : Ending position of match in database sequence
# 15 : ID