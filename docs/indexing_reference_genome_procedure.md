# INDEXING REFRENCE GENOME PROCEDURE

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

# Indexing reference genome for aligning
bwa index \
-p prefix/refence_genome \
refence_genome.fna ; \

# Adding fai index to reference
samtools faidx \
-o reference_genome.fna.fai \
reference_genome.fna

exit 0
```

[Next Procedure](data_preparation_procedure.md)
