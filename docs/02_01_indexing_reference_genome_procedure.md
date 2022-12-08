# INDEXING REFRENCE GENOME PROCEDURE

[Previous](02_00_procedure.md) | [Next](02_02_data_preparation_procedure.md)

## Quick Navigation

**[NOTEBOOK](../NOTEBOOK.md)**  
**[01 ECOGENETICS SETUP](01_00_ecogenetics_setup.md)**  
**[02 PROCEDURE](02_00_procedure.md)**  

- **[02 01 Indexing Reference Genome](02_01_indexing_reference_genome_procedure.md)**
- **[02 02 Data Preparation](02_02_data_preparation_procedure.md)**
- **[02 03 Initial Analysis Files](02_03_initial_analysis_procedure.md)**
- **[02 04 Genome Assembly and Annotation](02_04_genome_assembly_and_annotation.md)**

**[03 TERMINOLOGY](03_00_terminology.md)**  
**[04 SOFTWARE](04_00_software.md)**  
**[05 CLUSTER FUNCTIONS](05_00_cluster_functions.md)**

Needed software packages:

- bwa
- samtools

## Step 1: Indexing

### Script 01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 2:00:00

# Species specific reference genome, abosolute path (reference genome in FASTA format)
RG=""

# Indexing reference genome for aligning
bwa index \
-p "${RG%.*}" \
"$RG" ; \

# Adding fai index to reference
samtools faidx \
-o "$RG".fai \
"$RG"

exit 0
```

[Previous](02_00_procedure.md) | [Next](02_02_data_preparation_procedure.md)
