# INDEXING REFRENCE GENOME PROCEDURE

[Previous](02_00_procedure.md) | [Next](02_02_data_preparation_procedure.md)

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

Needed software packages:

- bwa
- samtools

## Step 0: Initialization

[Initialization script](../scripts/01_indexing_reference_genome/01_initialize_indexing_reference.sh)

## Step 1: Indexing

### Script 01

[01_01_indexing_reference.sh](../scripts/01_indexing_reference_genome/modules/01_01_indexing_reference.sh)  
Creates standard indexation files relating to BWA and a .fai index file from samtools.

[Previous](02_00_procedure.md) | [Next](02_02_data_preparation_procedure.md)
