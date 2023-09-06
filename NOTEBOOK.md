# **NOTEBOOK**

[Next](docs/01_00_ecogenetics_setup.md)

## Quick Navigation

**[NOTEBOOK](NOTEBOOK.md)**  
**[01 ECOGENETICS SETUP](docs/01_00_ecogenetics_setup.md)**  
**[02 PROCEDURE](docs/02_00_procedure.md)**  

- **[02 01 Indexing Reference Genome](docs/02_01_indexing_reference_genome_procedure.md)**
- **[02 02 Data Preparation](docs/02_02_data_preparation_procedure.md)**
- **[02 03 Initial Analysis Files](docs/02_03_initial_analysis_procedure.md)**
- **[02 04 Genome Assembly](docs/02_04_genome_assembly.md)**
- **[02 05 Genome Annotation](docs/02_05_genome_annotation.md)**

**[03 TERMINOLOGY](docs/03_00_terminology.md)**  
**[04 SOFTWARE](docs/04_00_software.md)**  
**[05 CLUSTER FUNCTIONS](docs/05_00_cluster_functions.md)**  

## **TO-DO**

### Might be useful

To polarize SNPs during the SFS estimation, we generated a genome consensus sequence for the okapi using ANGSD (flag -doFasta 1).

### For Museomics

Call VCF

```bash
freebayes  -f  reference.fa -L bamlist.txt> popname.vcf
```

filter to keep a minimum of 10x per genotype, no indels and no missing gentypes

```bash
vcftools --vcf AgUr_2022_J.vcf --recode --out name_DP10_noMISSING_noINDELS --minDP 10 --max-missing-count 0 --remove-indels
```

VCF -> SFS

```bash
conda install -c jaredgk -c bioconda py-popgen
```

Model-file

```bash
model_creator.py --model 1Pop --model-pop 1Pop name --pop-ind-file name name.txt
```

name.txt
  
- ind1  
- ind2  
- ind3  
- indn

BED-file  
GFF to BED

```bash
grep 'gene' gff_file > gff_file_gene

awk '{print $1, $4, $5}' gff_file_gene > gene.bed

bioawk  -v OFS='\t' -c fastx '{ print $name, length($seq) }' < sequences.fa > genome.txt

sort -k1,1 -k2,2n < gene.bed > gene_sorted.bed

sort -k1,1 -k2,2n < genome.txt > genome_sorted.txt

bedtools complement -i gene_sorted.bed -g genome_sorted.txt > intergenic.bed

vcf_to_sfs.py --vcf name_final_.recode.vcf --model-file out.model --modelname 1Pop  --folded --out name_sfs.txt
```

### Fra gff til intergene.bed som indeholder neutrale/intergene positioner

grep 'gene' gff_file > gff_file_gene

Tjek at det er kolonne 1,4 og 5 vi ønsker

```bash
awk '{print $1, $4, $5}' gff_file_gene > gene.bed
```

conda install -c bioconda bioawk

print the length of all scaffolds

```bash
bioawk  -v OFS='\t' -c fastx '{ print $name, length($seq) }' < genome_sequence.fa > genome.txt

sort -k1,1 -k2,2n < gene.bed > gene_sorted.bed

sort -k1,1 -k2,2n < genome.txt > genome_sorted.txt
```

conda install -c bioconda bedtools

```bash
bedtools complement -i gene_sorted.bed -g genome_sorted.txt > intergenic.bed
```

### Pileup -> vcf

```bash
bcftools call -vmO v -o <study.vcf>
```

---

After success of individual steps create pipelines to simplify workflow. For instance, combine mapping to reference genome and conversion to bam file, so there is no intermidiate storage of sam files

Always check for option to allocate multiple threads/cores to a job reduce computaion time

Look in to the use of the tee command for splitting output using pipes

### **Probably not going to use platypus**

    bwa index of reference
    samtools faidx of reference 

(used by platypus) index file need to be placed in the same directory as reference genome

    bwa mem -t -p | samtools sort -O BAM 

(samtools view step for file conversion not needed)

    samtools index
    platypus callVariants --bamFile=path/to/file.bam --refFile=path/to/referencefile.fa --output=path/to/output.vcf

creates vcf file from sorted .bam file and faidx indexed reference genome

### **Make mpileup file with bcftools**

    bcftools mpileup -Q 30 -q 30 -f path/to/referencegenome.fa > output.mpileup

`-Q` is cutoff value for baseQ
`-q` is cutoff value for mapQ

### **Site Frequency Spectrum (SFS)**

Possible how to:  
Use bcftools to create bcf file from sam then work on bcf file.
    bcftools mpileup -o out.bcf in.bam
    bcftools view -NIbl sites.file all.bcf > site.bcf
    bcftools view -cGP cond2 site.bcf > /dev/null 2> round1.afs
    bcftools view -cGP round1.afs site.bcf > /dev/null 2> round2.afs

Jinliang Yang Lab:  
ANGSD

    # write the bam files to a txt file
    mkdir bam_files
    mv sorted*.bam bam_files
    cd bam_files
    ls sorted*.bam > bam.txt

    samtools faidx Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa
    # run ANGSD to calculated folded SFS
    angsd -bam bamlist.txt -out output -doMajorMinor 1 -doMaf 1 -doSaf 2 -uniqueOnly 0 -anc Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa -minMapQ 30 -minQ 20 -nInd 20 -baq 1 -ref Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa -GL 1
    # use realSFS to calculate sfs
    realSFS output.saf.idx -fold 1 > output.sfs

In the local computer using R

    s <- scan('cache/output.sfs')
    s <- s[-c(1,length(s))]
    s <- s/sum(s)
    barplot(s,names=1:length(s), main='SFS')

### **X-forwarding**

Make X forwarding work to make it possible to use igv from the cluster

---

## **LOG**

### **31/08-2022**

BWA used to index Orchesella_cincta reference genome.  
Script: [01_index_reference_genome_orchesella_cincta.sh](scripts/01_index_reference_genome.sh)  
Output location: BACKUP/reference_genomes/orchesella_cincta/GCA_001718145.1

### **01/09-2022**

Mapping Orchesella_cincta sample to indexed reference genome using BWA  
First attempt surpassed reversed wall time of 4 hours and was thus cancelled prematurely. New run is has a time limit of 24 hours  
Script: [02_map_sample_to_reference.sh](scripts/02_map_sample_to_reference.sh)  
Standard BWA settings  
Output location: people/Jeppe_Bayer/steps/bwa/Orchesella_cincta/NYS-F  
Second attempt set to 24 hours also surpassed walltime

### **05/09-2022**

Mapping Orchesella_villosa? to Orchesella_cincta reference genome under the assumption that Orchesella_villosa and Orchesella_cincta are closely enough related.  
Script: [02_map_sample_to_reference.sh](scripts/02_map_sample_to_reference.sh)  
Standard BWA settings  
Output location: people/Jeppe_Bayer/steps/bwa/Orchesella_villosa?/FÅJ-C5

Conversion of BWA output file from .sam to .bam using samtool for Orchesella_cincta  
Script: [03_sam_to_bam.sh](scripts/03_sam_to_bam.sh)  
Flags used: -S specifying input as .sam, -b specifying output as .bam  
Output location: people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F

### **06/09-2022**

Previous attempts at mapping samples to reference genome were too resource restricted. Altered scripts resource demand to 48 hours, 6 cores, 6 G per core  
Follow-up on 07-09-2022

### **07/09-2022**

Mapping of both the Orchesella_cincta and Orchesella_villosa? sample to the Orchesella_cincta reference genome has been completed.  
Script: [02_map_sample_to_reference.sh](scripts/02_map_sample_to_reference.sh)  
Can be cut down to 2 GB per core. Mapping taking 18 to 21 hours producing sam files of 590 GB to 630 GB  
Output location: people/Jeppe_Bayer/steps/bwa/Orchesella_cincta/NYS-F and people/Jeppe_Bayer/steps/bwa/Orchesella_villosa?/FÅJ-C5 respectively

Converting both sam files to bam files. Estimated to take up much less space. After conversion the sam files will be deleted as not to use unnecessary space.  
Script: [03_sam_to_bam.sh](scripts/03_sam_to_bam.sh)  
Flags used: -b specifying output as .bam, -S has been deprecated.  
Output location: people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F and people/Jeppe_Bayer/steps/samtools/Orchesella_villosa?/FÅJ-C5 respectively  
Both file conversions took about 4 hours with single cores.

### **08/09-2022**

bam files have to be sorted as they are originally created from unsorted FASTQ files  
Sorting bam files, then trying to pipe them into coverage  
Checking coverage of mapped reads of Orchesella_cincta and Orchesella_villosa?  
Script: 04_sort_and_coverage.sh (deprecated)  
Sorted bam file is saved as `__.sorted.bam`  
In future pipe directly into sorting.  
Coverage output to unspecified format `.coverage`

### **09/09-2022**

Trying to gauge quality of maps and indexing .bam files with .bai files. .bai index files are used by other program such as IGV and other samtools commands to increase accessability and reduce time consumption  
.bai files are automatically accessed when need provided they are in the same directory as their partnered .bam file.  
Script: 06_flagstat.sh, 07_index_bam_to_bai.sh (deprecated, now [05_index_bam_file.sh](scripts/05_index_bam_file.sh) and [07_quality_check.sh](scripts/07_filtering_bam.sh))  
Used samtools idxstats to be able to get information on each contig, its length, number of mapped reads, and number of unmapped reads

Next up - Quality filtering, creating a new .bam file. When certain of the quality of new .bam file old .bam file can be deleted, as it can always be recomputed.  
When quality filtering is complete calculate allele frequency spectrum to look for species mixing in sample. AFS can be constructed on current data if quality filtering takes long.  
AFS might be possible to do using PoolHMM. PoolHMM needs data to be in mpileup file.  
Finish writing procedure section for flagstat, idxstats and coverage

Look for stop frameshift mutations - how?

### **12/09-2022**

Expanded NOTEBOOK with useful information about file formats, functions output and procedure  
Filtering BAM files for both Orchesella_cincta and Orchesella_villosa? to only include mapped, properly paired reads with MAPQ >=20  
Script: [07_filtering_bam.sh](scripts/07_filtering_bam.sh)

### **13/09-2022**

Counted number of reads in Orchesella_cincta BAM file after filtering and compared with unfiltered number of reads. Filtering left 78.8482 % of reads.  
Counted number of reads in Orchesella_villosa? BAM file after filtering and compared with unfiltered number of reads. Filtering left 27.2347 % of reads.  
Script: [08_check_reduction_in_reads.sh](scripts/08_check_reduction_in_reads.sh)

### **14/09-2022**

Cluster was down for maintanence. Looked into necessary data conversion for analysis. Platypus to create VCF files, bcftools to create mpileup.  
Created schematic of current workflow (see image under "Procedure").  
Extended information on files formats.  
Found several new useful articles

### **15/09-2022**

Changed NOTEBOOK to markdown format for better writing and reading experience.  
Created back-up of personal directory under the EcoGenetics project to GitHub <https://github.com/jeppebayer/EcoGenetics_Master>  
Needed to add platypus to conda environment. Incompatibility detected with version of python in environment ecogen. Created new environment based on compatibility with platypus. New environment named ecogen_prim and saved as [environmnet_primary.yml](environment_primary.yml). Old environment renamed to ecogen_secon and saved as [environment_secondary.yml](environment_secondary.yml)

### **16/09-2022**

Accidentally deleted notebook...  
Current version is recovered from GitHub and has lost whatever changes I've made today. Hopefully I think, THINK, that I only plotted out an optimized procedure... maybe... and a note on a job I'm trying to run

I want to redo some previous steps as I want to mark and remove duplicates. Trying to make a more efficient script that chains other scripts, e.g. each scripts end by queuing a new job containing the next script... I am having trouble. Maybe I should use gwf?

### **20/09-2022**

Script from 16/09 failed  
Rerun of [data preparations steps](scripts/data_preparation). Now fitted into fewer scripts which all get submitted at ones with dependencies so the start and of individual jobs is timed according to other jobs it depends on.  
samtools commands: faidx, collate, fixmate, markdup, have been incorporated into workflow.  
For markdup an optical duplicate distance of 100 has been chosen as the read length for sequencing was 150.  

### **21/09-2022**

Another failed run. Version of Samtools in use is 1.6 which does not support all the same functions and options as the most current version, several settings implemented in the script could not be use. Tried to alter jobs again, awaiting results.  
Presented outline for presentation for seminar day, got feedback.

### **22/09-2022**

Instructor work 10:00 - 14:00  
Working on presentation.

### **27/09-2022**

Still having trouble with marking of duplicates worked into pipe. Seems to be a problem in samtools 1.6? Have created a new environment without platypus, and thus newer version of everything else. Samtools is now version 1.15

### **28/09-2022**

Pipe still not working. Have rename script files with piped version to 00_brokenpipe1.sh and 00_brokenpipe2.sh. I have split the process in 11 seperate script jobs, one for each function/step. Running now. If it doesn't work might have to manually submit jobs and not do it in one batch with dependencies.

### **29/09-2022**

Problem has been identified and solved. Samtools fixmate function does NOT operate by writing to the standard output, thus any attempt at redirection or piping without filling the argument space used for output has led to empty named files.  
Deleted secondary environment, I don't think it will be relevant ever again. Also remove associated .yml file

### **04/10-2022**

No need to look for optical duplicates in markdup step, needed information is not available.  
Trying piped script on Ochesella_villosa?, mapping it to Orchesella_cincta reference, to test pipe and get an extra dataset to work with.  
Doing filtering of data. Removing duplicates and unmapped reads and retaining properly mapped reads with a MapQ $\ge$ 20.

### **05/10-2022**

Pipeline is working. Trying out data conversion. Working on making data availalbe for site frequency spectrum for sample data. Making variant call files. Ask Jesper about meaningful maxium read depth, highsest found ~19000

### **06/10-2022**

bcftools does not work with pooled sequences as bcftools only operates with ploidy of 1 or 2. Using samtools mpileup to make data available for popoolation, which cangive frequencies of pooled sequencing data. Popoolation should then be able to make it into a workable vcf file.

### **26/10-2022**

Continue work on 01_alt for adapterremoval. Ask Jesper about documentation for adapters used by bgi. Mads has experience using this program.

### **02/11 - 2022**

Over the last several days decisions have been made on standadizing procedures. I am now keeping an overview of the different procedures in markdown formatted document. These documents should include everything needed to replicated the scripts use for different procedures.

Standard directory construction for samples:

- BACKUP [directory]
  - population_genetics [directory]
    - <species_name> [directory]
      - <species_shorthand>_<location_shorthand> [directory]
        - sample_R1 [file]
        - sample_R2 [file]
        - analysis_ready_bam [file]

Standard directory construction for scripts:

- scripts [directory]
  - <procedure_number>_<procedure_name> [directory]
    - <procedure_number>_<script_number>_<procedure_name> (procedure name only on first script) [file]

### **07/11-2022**

When current running version of 02_data_preparation is done running:

- Change script directory to be built in and thus not need to be supplied.
- Change folder for final .bam file to be the population_genetics folder for the sample - same for the associated .bai file.
- Package all scripts relating to 02_data_preparation, with the exception of the init script, to be in a single folder called "core" or "modules" or something of the sort

### **17/11-2022**

Data preparation master script has been completed.  
Working on small overhaul to script 01 for indexing reference genomes.  
Working on script 03 to create initial analysis files for all population genetics samples. For every sample there is to be a pileup, VCF, SFS-complete, SFS-intergenic and SFS-nonsyn file. All are to be placed in the corresponding sample folder.  
Looking into the use of cactus for multiple sequence alignment. I have gotten a workflow from Jilong.  
Looking into genome assembly and annotation. I have gotten workflow files from Jilong.  

Need to do some general housekeeping and reorganization of files.  
Also have some pictures data need to made to schematics and need new illustrations for some proccesses.

### **18/11-2022**

Reduced datapreparation master script to a third of the size. Most likely also reducing computation time marginally.

### **24/11-2022**

Update on work:  
Data preparation master script has been completed.  
Working on small overhaul to script 01 for indexing reference genomes.  
Working on script 03 to create initial analysis files for all population genetics samples. For every sample there is to be a pileup, VCF, SFS-complete, SFS-intergenic and SFS-nonsyn file. All are to be placed in the corresponding sample folder.  
Looking into the use of cactus for multiple sequence alignment. I have gotten a workflow from Jilong.  
Looking into genome assembly and annotation. I have gotten workflow files from Jilong.  
--- For the museomics samples need to collect all samples within each species into one VCF file

Need to do some general housekeeping and reorganization of files.  
Also have some pictures data need to made to schematics and need new illustrations for some proccesses.

### **15/12-2022**

After correspondence with Simon Boitard, the creater of Pool-HMM, the issues relating to Pool-Hmm seem to possibly be related to the cluster itself. The program cannot currently be reliably used to create SFS for the samples but we still intend on using it for detecting possible selective sweeps. As an alternative to create SFS we, Jesper, Mads and I, have decided that I will use a small sub-program from Popoolation2 that creates files of the format sync. This is basically just a program to directly count variants occuring in a mpileup. I will then create a collection of scripts that use this format to create SFS files that can be used with Pool-HMM.
The created SFS will be folded and adhere to the following criteria:

- Not contain any reads with more than 2 alternatve bases
- Not contain any reads with ambiguous bases or deletions.
- Only contain reads with at least 2 observations for the MA
- All reads have a coverage of >= 200x
- Percentages will be rounded to nearest integer. Less than 0.5 percent is not counted

---

## **JOBS**

- **08/09-2022**  
    6414524  
    Orchesella_cincta sort

    6422047  
    Orchesella_villosa? sort and coverage  

    6422048  
    Orchesella_cincta coverage  

- **09/09-2022**  
    6441827  
    Orchesella_cincta flagstat

    6441836  
    Orchesella_villosa? flagstat

    6441875  
    Orchesella_cincta index of bam

    6441878  
    Orchesella_villosa? index of bam

    6441917  
    Orchesella_cincta idxstats

    6441918  
    Orchesella_villosa? idxstats

- **12/09-2022**  
    6511003  
    Orchesella_cincta filtering to only include mapped, properly paired reads with MAPQ >=20

    6511007  
    Orchesella_villosa? filtering to only include mapped, properly paired reads with MAPQ >=20

- **13/09-2022**  
    6530329  
    Orchesella_cincta - checking reduction in number of reads after filtering of BAM

    6530326  
    Orchesella_villosa? - checking reduction in number of reads after filtering of BAM

- **16/09-2022**  
    6583727  
    Job chain for Orchesella_cincta contains all data preparation steps except indexing of reference genome

- **20/09-2022**  
    6653874  
    Alternative to previous job chain
  - 6653875  
        Part one  
  - 6653876  
        Part two  
  - 6653877  
        Part three  
  - 6653878  
        Part four

- **22/09-2022**  
    6901513
    Main
  - 6901755
        Part two
  - 6901757
        Part three

- **24/09-2022**  
    6980166  
    intermediate  
  - 6980167  
        part two  
  - 6980168  
        part three

    6980169  
    no intermediate  
    - 6980170  
        part five  
    - 6980171  
        part six

- **27/09-2022**  
    7198212  
    init  
  - 7198213  
        part one  
  - 7198214  
        part two  
  - 7198216  
        part three  
  - 7198215  
        part five  
  - 7198217  
        part six

- **28/09-2022**  
    7223343  
    init  
  - 7223344  
        two. Time: 14:04:22  
  - 7223345  
        three. Time: 06:38:32  
  - 7223346  
        four. Time: 01:52:48  
  - 7223347  
        five. Failed due to typing error

- **29/09-2022**  
    7246667  
    init  
  - 7246672  
        five. Time: 00:56:49  
  - 7246673  
        six. Time: 05:03:41  
  - 7246674  
        seven. Failed due to typing error  

- **30/09-2022**  
    7276243  
    init  
  - 7276244  
        seven. Time: 04:44:24  
  - 7276245  
        eight. Time: 00:11:09  
  - 7276246  
        nine. Time: 00:12:24  
  - 7276247  
        ten. Time: 00:00:02  
  - 7276248  
        eleven. Time: 01:42:05  

- **03/10-2022**  
    7386389  
    12.sh creates depth file  
    probably not gonna use

- **04/10-2022**  
    7430351  
    init
  - 7430352  
        twelve. Time: 00:37:27  
  - 7430353  
        thirteen. Time: 00:16:59  

    7430354  
    alt_init  
    - 7430356  
        alt_two. Time: 20:52:45  
    - 7430357  
        alt_three. Time : 00:56:01  
    - 7430358  
        alt_four. Time: 00:20:54  

- **05/10-2022**  
    7447460  
    fourteen. bai index for filtered bam. Time: 00:05:26  
    7447498  
    ~~fifteen. make mpileup from filtered bam. Time: 07:20:47~~  
    7461944  
    ~~sixteen. creates variant call files from filtered bam. Time: 14:10:17~~  
    7461947  
    ~~alt_five.  bai index for filtered bam and creates variant call files from filtered bam. Time: 04:05:59~~  

- **06/10-2022**  
    7583734  
    fifteen. Time: 14:06:53  
    7583745  
    alt_five. Time. 04:27:30  

- **12/10-2022**  
    7870514  
    02 02 Orchesella_cincta - Output only zeroes

- **13/10-2022**  
    8111014  
    02 02 Orchesella_villosa?

- **14/10-2022**  
    8207478  
    02 02 Orchesella_cincta  
    8207609  
    02 03 Orchesella_villosa?  
    8209518  
    01 06_alt Orchesella_villosa?  
    8209533  
    01 06_alt Orchesella_cincta  

- **26/10-2022**  
    8517308  
    01_alt orchesella cincta

- **02/11-2022**
  8940003  
  01_01_indexing_reference_genome orchesella_cincta  
  8940740
  01_01_indexing_reference_genome pogonognathellus_flavescens

[Next](docs/01_00_ecogenetics_setup.md)