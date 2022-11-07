# **PROCEDURE**

## Procedure Steps

[01 Indexing Reference Genome](indexing_reference_genome_procedure.md)

[02 Data Preparation](data_preparation_procedure.md)

## Overview

(Needs to be overhauled)

1. Index reference genome using "bwa index"
2. Align sample sequences to reference genome using "bwa mem"
3. Convert SAM file created by "bwa mem" to BAM file using "samtools view"
4. Sort BAM file using "samtools sort"
5. Create index file for BAM file using "samtools index"
6. Quality check
   1. Create flagstat file for BAM file using "samtools falgstat"
   2. Create idxstats file for BAM file using "samtools idxstats"
   3. Create file of alignment coverage for BAM file using "samtools coverage"
7. Filter BAM file for low quality reads
8. Check number of remaining reads after filtering

![Workflow](../resources/workflow02.png)
