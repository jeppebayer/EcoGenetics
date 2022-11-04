# INITIAL ANALYSES PROCEDURE

[Previous Procedure](data_preparation_procedure.md)

Needed software packages:

- samtools

## Step 1: Create pileup file

### Script 01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

# Make mpileup file
samtools mpileup \
-C 50 \
-o sample.mpileup \
-f <reference_genome> \
filtered.bam

exit 0
```
