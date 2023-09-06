# **TERMINOLOGY**

[Previous](02_05_genome_annotation.md) | [Next](04_00_software.md)

## Quick Navigation

**[NOTEBOOK](../NOTEBOOK.md)**  
**[01 ECOGENETICS SETUP](01_00_ecogenetics_setup.md)**  
**[02 PROCEDURE](02_00_procedure.md)**  

- **[02 01 Indexing Reference Genome](02_01_indexing_reference_genome_procedure.md)**
- **[02 02 Data Preparation](02_02_data_preparation_procedure.md)**
- **[02 03 Initial Analysis Files](02_03_initial_analysis_procedure.md)**
- **[02 04 Genome Assembly](02_04_genome_assembly.md)**
- **[02 05 Genome Annotation](02_05_genome_annotation.md)**

**[03 TERMINOLOGY](03_00_terminology.md)**  
**[04 SOFTWARE](04_00_software.md)**  
**[05 CLUSTER FUNCTIONS](05_00_cluster_functions.md)**

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

[Previous](02_05_genome_annotation.md) | [Next](04_00_software.md)