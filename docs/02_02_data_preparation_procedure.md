# DATA PREPARATION PROCEDURE

[Previous](02_01_indexing_reference_genome_procedure.md) | [Next](02_03_initial_analysis_procedure.md)

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

- AdapterRemoval
- bwa
- samtools
- qualimap

Conda environmnet used:  
[environmnet_primary_from_history.yml](../environment_primary_from_history.yml)

## Step 0: Initialization

(Dependency tree slightly outdated)
![dependency_tree](../resources/dependency_chart_data_preparation.png)

## Step 1: Removal of overlapping sequences

### Script 01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

R1=

for R2 in "$sample"/*.fq.gz ; do
    if [ "$R1" ] ; then
        # 2nd entry
        AdapterRemoval \
        --threads 8 \
        --file1 "$R1" \
        --file2 "$R2" \
        --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        --minquality 25 \
        --minlength 20 \
        --basename "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed \
        --trimns \
        --trimqualities \
        --collapse

    else
        # First entry
        R1=$R2
    fi
done

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.discarded \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.settings \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.singleton.truncated

exit 0
```

### Script 02_01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Align sample to reference genome
bwa "$algo" -t "$cpus" \
"${RG%.*}" \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated \
| \

# Sort with regards to QNAME and convert to bam format
samtools sort -@ "$(("$cpus" - 1))" -n -O BAM \
-T "$WD"/temp/ \
-o "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
-

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair1.truncated \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.pair2.truncated

exit 0
```

### Script 02_02

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Concatenating collapsed single-end files
cat \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.collapsed \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.collapsed.truncated \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.all_collapsed

# Align sample to reference genome
bwa "$algo" -t "$cpus" \
"${RG%.*}" \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.all_collapsed \
| \

# Sort with regards to QNAME and convert to bam format
samtools sort -@ "$(("$cpus" - 1))" -n -O BAM \
-T "$WD"/temp/  \
-o "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_collapsed_aligned.bam \
-

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.collapsed \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.collapsed.truncated \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed.all_collapsed

exit 0
```

### Script 03

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Merge R1R2 and collapsed file to one name-sorted bam file
samtools merge -@ "$(("$cpus" - 1))" \
-o "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
-c -p -n \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_collapsed_aligned.bam

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_paired_aligned.bam \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_collapsed_aligned.bam

exit 0
```

## Step 2: Mark duplicates

### Script 04

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Add fixmate tag to alignment. Can only be done on name sorted alignment
samtools fixmate -@ "$(("$cpus" - 1))" -m -O BAM \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam \
- | \

# Position sort alignment and pipe output
samtools sort -@ "$(("$cpus" - 1))" -O BAM \
-T "$WD"/temp/ \
- | \

# Mark duplicates and save output as .coordinatesort.bam. Also outputs some realted statistics
samtools markdup -@ "$(("$cpus" - 1))" -s \
-f "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.markdupstats \
-T "$WD"/temp/ \
- \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_trimmed_complete_aligned.bam

exit 0
```

## Step 3: Statistics pre-filtering

### Script 05_01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates bai index for alignment
samtools index -@ "$(("$cpus" - 1))" -b \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam.bai ; \

# Creates idxstats file for alignment
samtools idxstats -@ "$(("$cpus" - 1))" \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.idxstats

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam.bai

exit 0
```

### Script 05_02

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates flagstat file for alignment 
samtools flagstat -@ "$(("$cpus" - 1))" \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.flagstat

exit 0
```

### Script 05_03

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates coverage file for alignment
samtools coverage \
-o "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.coverage \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

exit 0
```

### Script 05_04

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates stat file for alignment
samtools stats -@ "$(("$cpus" - 1))" \
-c 1,1000,1 \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/pre_filter_stats/"$(basename "$sample")"_markdup.stats

exit 0
```

## Step 4: Remove duplicates, unmapped reads and low quality mappings

### Script 06

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Removes duplicates and unmapped reads and keeps properly aligned reads with a MapQ >= 20
samtools view -b -@ "$(("$cpus" - 1))" \
-F 3844 -q 20 \
-o "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam && \

# Creates bai index for filtered alignment
samtools index -@ "$(("$cpus" - 1))" -b \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam.bai

exit 0
```

## Step 5: Statistics post-filtering

### Script 07_01

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates flagstat file for alignment 
samtools flagstat -@ "$(("$cpus" - 1))" \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.flagstat

exit 0
```

### Script 07_02

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates idxstats file for alignment
samtools idxstats -@ "$(("$cpus" - 1))" \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.idxstats

exit 0
```

### Script 07_03

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates coverage file for alignment
samtools coverage \
-o "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.coverage \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam

exit 0
```

### Script 07_04

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates text file listing number of reads in markdup file on first line
# number of reads in filtered file on second line and % remaining reads on the third line
samtools view -@ "$(("$cpus" - 1))" -c \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange && \
samtools view -@ "$(("$cpus" - 1))" -c \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
>> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange && \
awk \
'BEGIN{RS = "" ; FS = "\n"}{print "\n", $2/$1*100, "%"}' \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange \
>> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_markdup_to_filtered.readchange

# File removal
rm -f \
"$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/"$(basename "$sample")"_markdup.bam

exit 0
```

### Script 07_05

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

# Creates coverage file for alignment
samtools stats -@ "$(("$cpus" - 1))" \
"$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
> "$WD"/"$(basename "$script_path")"/"$(basename "$SD")"/"$(basename "$sample")"/post_filter_stats/"$(basename "$sample")"_filtered.stats

exit 0
```

### Script 07_06

```bash
#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

cpus=$1 # Number of CPUs
RG=$2 # Reference genome
SD=$3 # Species directory
WD=$4 # Working directory
sample=$5 # Sample directory
script_path=$6 # Path to script location
algo=$7 # Chosen algorithm

export DISPLAY=:0

# First checks whether a .gff file for the reference genome is available
for file in "$(dirname "$RG")"/*.gff; do
    if [ -e "$file" ]; then
        
        # .gff file is available
        gff=$file

        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -gff "$gff" \
        -outdir "$SD"/"$(basename "$sample")"/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=20G
        exit 0

    else

        # .gff file is not available
        qualimap bamqc \
        -bam "$SD"/"$(basename "$sample")"/"$(basename "$sample")"_filtered.bam \
        -outdir "$SD"/"$(basename "$sample")"/qualimap \
        -outfile "$(basename "$sample")"_qualimap.pdf \
        -outformat PDF \
        --java-mem-size=20G
        exit 0
    fi
done
```

[Previous](02_01_indexing_reference_genome_procedure.md) | [Next](02_03_initial_analysis_procedure.md)
