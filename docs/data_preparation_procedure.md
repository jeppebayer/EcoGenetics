# DATA PREPARATION PROCEDURE

[Previous Procedure](indexing_reference_genome_procedure.md)

Needed software packages:

- AdapterRemoval
- bwa
- samtools
- awk

## Step 0: Initialization

![dependency_tree](../resources/dependency_chart_data_preparation.png)

### Script 00

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH --cpus-per-task 1
#SBATCH --time 00:10:00

# ----------------- Configuration ----------------------------------------

# Species specific reference Genome Directory
RGD=""

# Species specific sample Directory
SD=""

# Working Directory
WD=""

# ----------------- Script Queue -----------------------------------------

# AdapterRemoval
jid1=$(sbatch --parsable 01.sh)

# Aligning to reference
jid2_1=$(sbatch --parsable --dependency=afterok:"$jid1" 02_1.sh)
jid2_2=$(sbatch --parsable --dependency=afterok:"$jid1" 02_2.sh)

# Merging of alignment files
jid3=$(sbatch --parsable --dependency=afterok:"$jid2_1":"$jid2_2" 03.sh)

# Marking duplicates
jid4=$(sbatch --parsable --dependency=afterok:"$jid3" 04.sh)

# Statistics pre-filtering
jid5_1=$(sbatch --parsable --dependency=afterok:"$jid4" 05_1.sh)
jid5_2=$(sbatch --parsable --dependency=afterok:"$jid4" 05_2.sh)
jid5_3=$(sbatch --parsable --dependency=afterok:"$jid4" 05_3.sh)
jid5_4=$(sbatch --parsable --dependency=afterok:"$jid4" 05_4.sh)

# Removal of duplicates, unmapped reads and low quality mappings
jid6=$(sbatch --parsable --dependency=afterok:"$jid5_1":"$jid5_2":"$jid5_3":"$jid5_4" 06.sh)

# Statistics post-filtering
jid7_1=$(sbatch --parsable --dependency=afterok:"$jid6" 07_1.sh)
jid7_2=$(sbatch --parsable --dependency=afterok:"$jid6" 07_2.sh)
jid7_3=$(sbatch --parsable --dependency=afterok:"$jid6" 07_3.sh)
jid7_4=$(sbatch --parsable --dependency=afterok:"$jid6" 07_4.sh)
jid7_5=$(sbatch --parsable --dependency=afterok:"$jid6" 07_5.sh)

exit 0
```

## Step 1: Removal of overlapping sequences

### Script 01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 12:00:00

AdapterRemoval \
--threads 8 \
--file1 <R1> \
--file2 <R1> \
--adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
--adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
--minquality 25 \
--minlength 20 \
--basename trimmed \
--trimns \
--trimqualities \
--collapse

exit 0
```

### Script 02_01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 30:00:00

# Align sample to reference genome
bwa mem -t 8 \
<reference_genome> \
<trimmed_R1> <trimmed_R2> \
| \

# Sort with regards to QNAME and convert to bam format
samtools sort -@ 7 -n -O BAM \
-T path/to/temp_folder/ \
-o trimmed_R1R2_aligned.bam \
-

exit 0
```

### Script 02_02

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 30:00:00

cat trimmed_collapsed trimmed_collapsed_truncated > all_collapsed

# Align sample to reference genome
bwa mem -t 8 \
<reference_genome> \
<all_collapsed> \
| \

# Sort with regards to QNAME and convert to bam format
samtools sort -@ 7 -n -O BAM \
-T path/to/temp_folder/ \
-o trimmed_collapsed_aligned.bam \
-

exit 0
```

### Script 03

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 36:00:00

# Merge R1R2 and collapsed file to one name-sorted bam file
samtools merge -@ 7 \
-o complete_sample_aligned.bam \
-c -p -n \
trimmed_R1R2_aligned.bam \
trimmed_collapsed_aligned.bam

exit 0
```

## Step 2: Mark duplicates

### Script 04

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 30:00:00

# Add fixmate tag to alignment. Can only be done on name sorted alignment
samtools fixmate -@ 7 -m -O BAM \
complete_sample_aligned.bam \
- | \

# Position sort alignment and pipe output
samtools sort -@ 7 -O BAM \
-T path/to/temp_folder/ \
- | \

# Mark duplicates and save output as .coordinatesort.bam. Also outputs some realted statistics
samtools markdup -@ 7 -s \
-f markdup.markdupstats \
-T path/to/temp_folder/ \
- \
markdup.bam

exit 0
```

## Step 3: Statistics pre-filtering

### Script 05_01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Creates bai index for alignment
samtools index -@ 7 -b \
.markdup.bam \
> markdup.bam.bai ; \

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
markdup.bam \
> markdup.idxstats

exit 0
```

### Script 05_02

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Creates flagstat file for alignment 
samtools flagstat -@ 7 \
markdup.bam \
> markdup.flagstat

exit 0
```

### Script 05_03

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Creates coverage file for alignment
samtools coverage \
-o markdup.coverage \
markdup.bam

exit 0
```

### Script 05_04

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Creates coverage file for alignment
samtools stats -@ 7 \
-c 1,1000,1 \
markdup.bam \
> markdup.stats

exit 0
```

## Step 4: Remove duplicates, unmapped reads and low quality mappings

### Script 06

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 24:00:00

# Removes duplicates, unmapped reads, supplementary alignments, secondary alignments and  reads with a MapQ >= 20
samtools view -b -@ 7 \
-F 3844 -q 20 \
-o filtered.bam \
markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ 7 -b \
filtered.bam \
> filtered.bam.bai

exit 0
```

## Step 5: Statistics post-filtering

### Script 07_01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 01:00:00

# Creates flagstat file for alignment 
samtools flagstat -@ 7 \
filtered.bam \
> filtered.flagstat

exit 0
```

### Script 07_02

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 01:00:00

# Creates idxstats file for alignment
samtools idxstats -@ 7 \
filtered.bam \
> filtered.idxstats

exit 0
```

### Script 07_03

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 01:00:00

# Creates coverage file for alignment
samtools coverage \
-o filtered.coverage \
filtered.bam

exit 0
```

### Script 07_04

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 06:00:00

# Creates text file listing number of reads in markdup file on first line
# number of reads in filtered file on second line and % remaining reads on the third line
samtools view -@ 7 -c \
markdup.bam \
> markdup_to_filtered.readchange && \
samtools view -@ 7 -c \
filtered.bam \
>> markdup_to_filtered.readchange && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
markdup_to_filtered.readchange \
>> markdup_to_filtered.readchange

exit 0
```

### Script 07_05

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 1:00:00

# Creates coverage file for alignment
samtools stats -@ 7 \
-c 1,1000,1 \
filtered.bam \
> filtered.stats

exit 0
```

[Next Procedure](initial_analysis_procedure.md)
