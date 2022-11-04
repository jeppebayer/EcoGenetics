# **TERMINOLOGY**

## **MapQ**

Mapping quality. A measure of how likely a given sequence alignment is to be located were its reported location is.  
If a read's mapping quality is low, especially if it's 0, the read maps to multiple locations on the genome (they are multi-hit or multi-mapping reads), and one cannot be sure wether the reported location is the correct one.  
https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/

## **$\pi_N/\pi_S$**

The ratio of non-synonymous to synonymous site diversity.

## **Distribution of Fitness Effects (DFE)**

Describes what proportion of new mutations are advantageous, neutral or deleterious.

## **Genomc Evolutionary Rate Profiling (GERP)**

Identifies constrained elements in a multiple sequence alignment by quantifying substitution deficits (substitutions that would have occurred if the element were neutral). Rejected substitutions are a natural measure of constraint that reflects the strength of past purifying selection.

## **File formats**

**.SAM**  
format for storing read alignments  
row-oriented tab-delimited text file

**.BAM**  
binary counterpart to .SAM. contains same information, consumes much less space. Allows more efficient processing of data

**.VCF**  
format for storing genetic variants  
row-oriented tab-delimited text file

**.BCF**  
binary counterpart to .VCF. contains same information, consumes much less space. Allows more efficient processing of data