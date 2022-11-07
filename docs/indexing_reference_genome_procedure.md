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

[Next Procedure](data_preparation_procedure.md)
