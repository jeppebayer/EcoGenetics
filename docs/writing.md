# **Overview**

- [**Overview**](#overview)
  - [**Methods**](#methods)
    - [**Mapping and Filtering**](#mapping-and-filtering)
      - [***Faster X***](#faster-x)
      - [***PopGen***](#popgen)
      - [***Museomics***](#museomics)
  - [**Acknowledgements**](#acknowledgements)
  - [**References**](#references)
    - [Faster X](#faster-x-1)
    - [AdapterRemoval](#adapterremoval)
    - [BWA](#bwa)
    - [BWA settings for historic samples](#bwa-settings-for-historic-samples)
    - [SAMtools](#samtools)
    - [Qualimap 2](#qualimap-2)
    - [HiFiAdapterFilt](#hifiadapterfilt)
    - [Jellyfish](#jellyfish)
    - [GenomeScope2](#genomescope2)
    - [hifiasm](#hifiasm)
    - [BUSCO](#busco)
    - [purge\_dups](#purge_dups)
    - [minimap2](#minimap2)
    - [DIAMOND](#diamond)
    - [BLAST](#blast)
    - [NCBI databases](#ncbi-databases)
    - [Uniprot database](#uniprot-database)
    - [BlobTools](#blobtools)
    - [RepeatModeler2](#repeatmodeler2)
    - [RepeatMasker](#repeatmasker)
    - [GIRI RepBase database](#giri-repbase-database)
    - [BEDTools](#bedtools)
    - [Freebayes](#freebayes)
    - [SnpEff](#snpeff)
    - [SnpSift](#snpsift)
    - [AGAT](#agat)
    - [Progressive Cactus](#progressive-cactus)
    - [HAL](#hal)
  - [**Code**](#code)
    - [**Utilities**](#utilities)
    - [**01 - Indexation of Reference Genome**](#01---indexation-of-reference-genome)
      - [**Conda Enviroment**](#conda-enviroment)
      - [**Indexing reference genome**](#indexing-reference-genome)
    - [**02 - Data Preparation, Mapping and Filtering**](#02---data-preparation-mapping-and-filtering)
      - [**Conda Enviroment**](#conda-enviroment-1)
      - [**AdapterRemoval**](#adapterremoval-1)
      - [**Alignment**](#alignment)
      - [**----- Paired-end (Faster X, PopGen, Modern)**](#------paired-end-faster-x-popgen-modern)
      - [**---------- Non-collapsed**](#-----------non-collapsed)
      - [**---------- Collapsed**](#-----------collapsed)
      - [**--------------- Align sample to reference genome**](#----------------align-sample-to-reference-genome)
      - [**--------------- Merge non-collapsed and collapsed alignment into one name-sorted BAM file**](#----------------merge-non-collapsed-and-collapsed-alignment-into-one-name-sorted-bam-file)
      - [**----- Single-end (Historic)**](#------single-end-historic)
      - [**Mark duplicates**](#mark-duplicates)
      - [**Filtering**](#filtering)
      - [**Extract unmapped reads**](#extract-unmapped-reads)
      - [**Post filtering stats**](#post-filtering-stats)
    - [**03 - Initial Analysis Files, VCF and SFS**](#03---initial-analysis-files-vcf-and-sfs)
      - [**Conda Enviroment**](#conda-enviroment-2)
      - [**Create VCF**](#create-vcf)
      - [**Create SnpEff Database**](#create-snpeff-database)
      - [**Annotation with SnpEff**](#annotation-with-snpeff)
    - [**04 - Genome Assembly**](#04---genome-assembly)
      - [**Conda Enviroment**](#conda-enviroment-3)
      - [**HiFi adapter filtering**](#hifi-adapter-filtering)
      - [**K-mer analysis**](#k-mer-analysis)
      - [**HiFi assembly**](#hifi-assembly)
      - [**----- Making a FASTA formatted file from GFA**](#------making-a-fasta-formatted-file-from-gfa)
      - [**BUSCO analysis**](#busco-analysis)
      - [**Purge dups**](#purge-dups)
      - [**BlobTools**](#blobtools-1)
      - [**----- COV file**](#------cov-file)
      - [**----- HITS files**](#------hits-files)
      - [**----- Blobplot**](#------blobplot)
    - [**05 - Repeat Detection, Annotating and Masking**](#05---repeat-detection-annotating-and-masking)
      - [**Building database for RepeatModeler2**](#building-database-for-repeatmodeler2)
      - [**Creating model with RepeatModeler2**](#creating-model-with-repeatmodeler2)
      - [**Repeat masking based on GIRI RepBase**](#repeat-masking-based-on-giri-repbase)
      - [**Repeat masking based on RepeatModeler2 model**](#repeat-masking-based-on-repeatmodeler2-model)
      - [**Combining results from repeat masking**](#combining-results-from-repeat-masking)
      - [**Complete soft masking of genome assembly and create BED file of repeat regions**](#complete-soft-masking-of-genome-assembly-and-create-bed-file-of-repeat-regions)
    - [**Miscellaneous**](#miscellaneous)
      - [**Make DIAMOND database from UniProt Reference Proteomes**](#make-diamond-database-from-uniprot-reference-proteomes)
      - [**Setting up Progressive cactus aligner**](#setting-up-progressive-cactus-aligner)

## **Methods**

### **Mapping and Filtering**

Prior to mapping, the re-sequencing data processed with AdapterRemoval (Schubert et al. 2016) to search for and remove any residuals of known adapter sequences. Furthermore, reads were trimmed at the 5' and 3' termini for Ns and bases with a quality lower than 25, reads shorter than 20 base-pairs were also discarded, and for paired-end reads overlapping mates were merged into a single read and the base quality was recalculated if the overlap was at least 11 base pairs in length (default) with a maximum mismatch ratio of 1/3 (default).

#### ***Faster X***

Polymoirphisms - Mapping and filtering

Re-sequencing data was first processed with AdapterRemoval version 2.3.2 (Schubert et al. 2016) to search for and remove any residuals of known adapter sequences. Furthermore, reads were trimmed at the 5’ and 3’ termini for ambiguous bases and bases with a quality lower than 25, reads shorter than 20 base-pairs were also discarded, and for overlapping mates were merged into a single read with a recalculated base quality if the overlap was at least 11 base-pairs in length with a maximum mismatch ratio of 1/3. 

Re-sequencing data from each species was mapped to their reference genomes using the BWA-MEM algorithm in BWA version 0.7.17-r1188 (Li 2013). The alignments were then converted to BAM format and name sorted using SAMtools version 1.16.1 (Li et al. 2009; Danecek et al. 2021). To identify duplicates the alignments were processed with ‘samtools fixmate’ and then ‘samtools markdup’. Each alignment was filtered using ‘samtools view’ removing alignments flagged as any of the following: unmapped reads, not primary alignment, failed reads, optical or PCR duplicates, supplementary reads, or alignments with a mapping quality of less than 20. Alignments were evaluated using Qualimap 2 version 2.2.2-dev (Okonechnikov et al. 2016).

#### ***PopGen***

Each pooled re-sequencing data sample file was mapped to a reference genome of their respective species using the BWA-MEM algorithm from BWA (Li H. 2013).

#### ***Museomics***

For modern samples (2022 or newer) the re-sequencing data was mapped to reference genomes using the BWA-MEM algorithm from BWA (Li H. 2013).
For historic samples (older than 2022) the re-sequencing data was mapped to reference genomes using the BWA-backtrack algorithm from BWA (Li H. and Durbin R. 2009).

## **Acknowledgements**

All of the computing for this project was performed on the GenomeDK cluster. We would like to thank GenomeDK and Aarhus University for providing computational resources and support that contributed to these research results.

## **References**

### Faster X

- Jesper Bechsgaard and others, Evidence for Faster X Chromosome Evolution in Spiders, *Molecular Biology and Evolution*, Volume 36, Issue 6, June 2019, Pages 1281–1293, <https://doi.org/10.1093/molbev/msz074>
- Maria Izabel A Cavassim and others, Recombination Facilitates Adaptive Evolution in Rhizobial Soil Bacteria, *Molecular Biology and Evolution*, Volume 38, Issue 12, December 2021, Pages 5480–5490, <https://doi.org/10.1093/molbev/msab247>
- Al-saffar & Hahn 2022, Evaluating methods for estimating the proportion of adaptive amino acid substitutions, <https://www.biorxiv.org/content/10.1101/2022.08.15.504017v1.full.pdf>
- Judith E. Mank and others, EFFECTIVE POPULATION SIZE AND THE FASTER-X EFFECT: EMPIRICAL RESULTS AND THEIR INTERPRETATION, *Evolution*,  Volume 64, Issue 3, 1 March 2010, Pages 663–674, <https://doi.org/10.1111/j.1558-5646.2009.00853.x>
- Beatriz Vicoso , Brian Charlesworth, EFFECTIVE POPULATION SIZE AND THE FASTER-X EFFECT: AN EXTENDED MODEL, *Evolution*, Volume 63, Issue 9, 1 September 2009, Pages 2413–2426, <https://doi.org/10.1111/j.1558-5646.2009.00719.x>

### AdapterRemoval

- Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. *BMC Research Notes*, 12;9(1):88. <http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2>

### BWA

- Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*, 25, 1754-1760. [PMID: 19451168]. (if you use the BWA-backtrack algorithm)
- Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN]. (if you use the BWA-MEM algorithm or the fastmap command, or want to cite the whole BWA package)

### BWA settings for historic samples

- Palkopoulou E, Mallick S, Skoglund P, Enk J, Rohland N, Li H, Omrak A, Vartanyan S, Poinar H, Götherström A, Reich D, Dalén L. Complete genomes reveal signatures of demographic and genetic declines in the woolly mammoth. *Curr Biol.* 2015 May 18;25(10):1395-400. doi: 10.1016/j.cub.2015.04.007. Epub 2015 Apr 23. PMID: 25913407; PMCID: PMC4439331.

### SAMtools

(BCFtools)

- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, Twelve years of SAMtools and BCFtools, *GigaScience* (2021) 10(2) giab008 [33590861]
- Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, *Bioinformatics* (2009) 25(16) 2078-9 [19505943]

### Qualimap 2

- Konstantin Okonechnikov and others, Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data, *Bioinformatics*, Volume 32, Issue 2, January 2016, Pages 292–294, <https://doi.org/10.1093/bioinformatics/btv566>

### HiFiAdapterFilt

- Sim, Sheina B. and Corpuz, Renee L. and Simmonds, Tyler J. and Geib, Scott M., (2022) HiFiAdapterFilt, a memory efficient read processing pipeline, prevents occurrence of adapter sequence in PacBio HiFi reads and their negative impacts on genome assembly

### Jellyfish

- Guillaume Marcais and Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. *Bioinformatics* (2011) 27(6): 764-770 (first published online January 7, 2011) doi:10.1093/bioinformatics/btr011

### GenomeScope2

- Ranallo-Benavidez, T.R., Jaron, K.S. & Schatz, M.C. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* 11, 1432 (2020). <https://doi.org/10.1038/s41467-020-14998-3>
- Vurture, GW, Sedlazeck, FJ, Nattestad, M, Underwood, CJ, Fang, H, Gurtowski, J, Schatz, MC (2017) *Bioinformatics* doi: <https://doi.org/10.1093/bioinformatics/btx153>

### hifiasm

- Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., Li H. (2021) Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nat Methods*, 18:170-175. <https://doi.org/10.1038/s41592-020-01056-5>
- Cheng, H., Jarvis, E.D., Fedrigo, O., Koepfli, K.P., Urban, L., Gemmell, N.J., Li, H. (2022) Haplotype-resolved assembly of diploid genomes without parental data. *Nature Biotechnology*, 40:1332–1335. <https://doi.org/10.1038/s41587-022-01261-x>

### BUSCO

- Manni M., Berkeley M.R., Seppey M., Simao F.A., Zdobnov E.M. 2021. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. arXiv:2106.11799 [q-bio] [Internet]. Available from: <http://arxiv.org/abs/2106.11799>

### purge_dups

- Guan D, McCarthy SA, Wood J, Howe K, Wang Y, Durbin R. Identifying and removing haplotypic duplication in primary genome assemblies. *Bioinformatics*. 2020 May 1;36(9):2896-2898. doi: 10.1093/bioinformatics/btaa025. PMID: 31971576; PMCID: PMC7203741.

### minimap2

- Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*. 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191. PMID: 29750242; PMCID: PMC6137996.

### DIAMOND

- Buchfink, B., Reuter, K. & Drost, HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods* 18, 366–368 (2021). <https://doi.org/10.1038/s41592-021-01101-x>

### BLAST

- Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. *J Mol Biol.* 1990 Oct 5;215(3):403-10. doi: 10.1016/S0022-2836(05)80360-2. PMID: 2231712.
- etc. <https://blast.ncbi.nlm.nih.gov/doc/blast-help/references.html#references>

### NCBI databases

- Sayers EW, Bolton EE, Brister JR, Canese K, Chan J, Comeau DC, Connor R, Funk K, Kelly C, Kim S, Madej T, Marchler-Bauer A, Lanczycki C, Lathrop S, Lu Z, Thibaud-Nissen F, Murphy T, Phan L, Skripchenko Y, Tse T, Wang J, Williams R, Trawick BW, Pruitt KD, Sherry ST. Database resources of the national center for biotechnology information. *Nucleic Acids Res.* 2022 Jan 7;50(D1):D20-D26. doi: 10.1093/nar/gkab1112. PMID: 34850941; PMCID: PMC8728269.

### Uniprot database

- The UniProt Consortium, UniProt: the Universal Protein Knowledgebase in 2023, *Nucleic Acids Research*, Volume 51, Issue D1, 6 January 2023, Pages D523–D531, <https://doi.org/10.1093/nar/gkac1052>

### BlobTools

- Laetsch DR and Blaxter ML. BlobTools: Interrogation of genome assemblies [version 1; peer review: 2 approved with reservations]. F1000Research 2017, 6:1287 (<https://doi.org/10.12688/f1000research.12232.1>)

### RepeatModeler2

- Flynn JM, Hubley R, Goubert C, Rosen J, Clark AG, Feschotte C, Smit AF. RepeatModeler2 for automated genomic discovery of transposable element families. *Proc Natl Acad Sci U S A*. 2020 Apr 28;117(17):9451-9457. doi: 10.1073/pnas.1921046117. Epub 2020 Apr 16. PMID: 32300014; PMCID: PMC7196820.

### RepeatMasker

- Tarailo-Graovac M, Chen N. Using RepeatMasker to identify repetitive elements in genomic sequences. *Curr Protoc Bioinformatics*. 2009 Mar;Chapter 4:4.10.1-4.10.14. doi: 10.1002/0471250953.bi0410s25. PMID: 19274634.

### GIRI RepBase database

- Bao, W., Kojima, K.K., Kohany, O. Repbase Update, a database of repetitive elements in eukaryotic genomes. *Mob DNA*, 2015;6:11

### BEDTools

- Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Jan 28. PMID: 20110278; PMCID: PMC2832824.

### Freebayes

- Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. *arXiv preprint arXiv:1207.3907 [q-bio.GN]* 2012

### SnpEff

- "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672

### SnpSift

- "Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift", Cingolani, P., et. al., Frontiers in Genetics, 3, 2012.

### AGAT

- Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. <https://www.doi.org/10.5281/zenodo.3552717>

### Progressive Cactus

- Armstrong, J., Hickey, G., Diekhans, M. et al. Progressive Cactus is a multiple-genome aligner for the thousand-genome era. *Nature* 587, 246–251 (2020). <https://doi.org/10.1038/s41586-020-2871-y>

### HAL

- Glenn Hickey and others, HAL: a hierarchical format for storing and analyzing multiple genome alignments, *Bioinformatics*, Volume 29, Issue 10, May 2013, Pages 1341–1342, <https://doi.org/10.1093/bioinformatics/btt128>

## **Code**

-----

### **Utilities**

- Ripgrep
- seqkit
- bioawk

### **01 - Indexation of Reference Genome**

-----

#### **Conda Enviroment**

Date: 01/06-2023

```yaml
name: data_prep
channels:
  - genomedk
  - gwforg
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - adapterremoval=2.3.2
  - qualimap=2.2.2d
  - matplotlib
  - ripgrep
  - bioawk
  - samtools=1.16.1
  - perl
  - pandas
  - bwa=0.7.17
prefix: /home/jepe/miniconda3/envs/data_prep
```

#### **Indexing reference genome**

```bash
# Index database sequences in the FASTA format
bwa index \
    # Prefix of the output database. Default 'same as db filename'
    -p <reference_genome_prefix> \
    <reference_genome>

# Index reference sequence in the FASTA format.
samtools faidx \
    # Write to file rather than to stdout.
    -o <reference_genome>.fai \
    <reference_genome>
```

### **02 - Data Preparation, Mapping and Filtering**

-----

#### **Conda Enviroment**

Date: 01/06-2023

```yaml
name: data_prep
channels:
  - genomedk
  - gwforg
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - adapterremoval=2.3.2
  - qualimap=2.2.2d
  - matplotlib
  - ripgrep
  - bioawk
  - samtools=1.16.1
  - perl
  - pandas
  - bwa=0.7.17
prefix: /home/jepe/miniconda3/envs/data_prep
```

#### **AdapterRemoval**

```bash
# (can use '--threads' for multithreading)
AdapterRemoval \
    --file1 <R1> \
    # Not present when single-end
    --file2 <R2> \
    # Adapters used by BGI
    --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
    # Not present when single-end. Adapter used by BGI
    --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
    # Set the threshold for trimming low quality bases using --trimqualities and --trimwindows. Default is 2.
    --minquality 25 \
    # Reads shorter than this length are discarded following trimming. Defaults to 15.
    --minlength 20 \
    --basename <basename> \
    # Trim consecutive Ns from the 5' and 3' termini.  If  quality  trimming  is  also  enabled  (--trimqualities),  then stretches of mixed low-quality bases and/or Ns are trimmed.
    --trimns \
    # Trim  consecutive  stretches  of  low  quality bases (threshold set by --minquality) from the 5' and 3' termini. If trimming of Ns is also enabled (--trimns), then stretches of mixed low-quality bases and Ns are trimmed.
    --trimqualities \
    # Not present when single-end. In paired-end mode, merge overlapping mates into a single and recalculate the quality scores. In  single-end  mode, attempt to identify templates for which the entire sequence is available. In both cases, complete "collapsed" reads are written with a 'M_' name prefix, and "collapsed" reads which are trimmed due to quality  settings  are  written with  a 'MT_' name prefix. The overlap needs to be at least --minalignmentlength (default 11) nucleotides, with a maximum number of mismatches determined by --mm (default maximum mismatch rate of 1/3).
    --collapse
```

<u>*Output files*</u>

For paired-end reads:

- File containing mate 1 reads  
  `<basename>.pair1.truncated`  
- File containing mate 2 reads  
  `<basename>.pair2.truncated`  
- File containing overlapping mate-pairs which have been merged into a single read  
  `<basename>.collapsed`  
- File containing collapsed reads which were trimmed due to the presence of low-quality or ambiguous nucleotides  
  `<basename>.collapsed.truncated`  
- File containing reads discarded due to the `--minlength` option  
  `<basename>.discarded`  
- File containing information on the parameters used in the run as well as overall statistics on the reads after trimming  
  `<basename>.settings`  
- File containing paired reads for which the mate has been discarded  
  `<basename>.singleton.truncated`  

For single-end reads:

- File containing reads  
  `<basename>.truncated`  
- File containing reads discarded due to the `--minlength` option  
  `<basename>.discarded`  
- File containing information on the parameters used in the run as well as overall statistics on the reads after trimming  
  `<basename>.settings`  

#### **Alignment**

#### **----- Paired-end (Faster X, PopGen, Modern)**

#### **---------- Non-collapsed**

```bash
# (can use '-t' for multithreading)
bwa mem \
    # Adds 'read-group' header line based on sample name
    -R "@RG\tID:<sample_name>\tSM:<sample_name>" \
    <reference_genome_prefix> \
    <basename.pair1.truncated> \
    <basename.pair2.truncated> \
# (can use '-@' for multithreading)
| samtools sort \
    # Sort by read name
    -n \
    # Output in BAM format
    -O BAM \
    -T <temporary_file_directory> \
    -o <paired_aligned.bam> \
    -
```

#### **---------- Collapsed**

#### **--------------- Align sample to reference genome**

```bash
# (can use '-t' for multithreading)
bwa mem \
    # Adds 'read-group' header line based on sample name
    -R "@RG\tID:<sample_name>\tSM:<sample_name>" \
    <reference_genome_prefix> \
    <(cat \
        <basename.collapsed> \
        <basename.collapsed.truncated>) \
# (can use '-@' for multithreading)
| samtools sort \
    # Sort by read name
    -n \
    # Output in BAM format
    -O BAM \
    -T <temporary_file_directory> \
    -o <collapsed_aligned.bam> \
    -
```

#### **--------------- Merge non-collapsed and collapsed alignment into one name-sorted BAM file**

```bash
# Merges multiple sorted alignment files, producing a single sorted output file that contains all the input records and maintains the existing order. (can use '-@' for multithreading)
samtools merge \
    # Keep only one @RG header line in the output
    -c \
    # Keep only one @PG header line in the output
    -p \
    # Indicates that input files are sorted by read names rather than chromosomal coordinates
    -n \
    -o <complete_aligned.bam> \
    <paired_aligned.bam> \
    <collapsed_aligned.bam>
```

#### **----- Single-end (Historic)**

Settings used with BWA, chosen from () as the had proved useful for ancient/historic DNA.

"Merged reads were mapped against the reference genome using parameters -l 16500 –n 0.01 –o 2, which deactivate seeding, allow more substitutions and permit up to two gaps (instead of one), using BWA’s ’aln’ algorithm to construct suffix arrays"

```bash
# Finds Suffix Array (SA) coordinates for input reads. (can use '-t' for multithreading)
bwa aln \
    # Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. Default inf
    -l 16500 \
    # Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. Default 0.04
    -n 0.01 \
    # Maximum number of gap opens. Default 1
    -o 2 \
    <reference_genome_prefix> \
    <basename.truncated> \
# Generate alignments in the SAM format given single-end reads. Repetitive hits will be randomly chosen.
| bwa samse \
    # Adds 'read-group' header line based on sample name
    -r "@RG\tID:<sample_name>\tSM:<sample_name>" \
    <reference_genome_prefix> \
    - \
    <basename.truncated> \
# (can use '-@' for multithreading)
| samtools sort \
    # Sort by read name
    -n \
    # Output in BAM format
    -O BAM \
    -T <temporary_file_directory> \
    -o <complete_aligned.bam> \
    -
```

#### **Mark duplicates**

```bash
# Fills in mate coordinates, ISIZE and mate related flags on a name sorted alignment. (can use '-@' for multithreading)
samtools fixmate \
    # Adds 'mate score' tag. Used by 'samtools markdup' to select the best reads to keep.
    -m \
    # Output in BAM format
    -O BAM \
    <complete_aligned.bam> \
    - \
# Position sort alignment by leftmost coordinate. (can use '-@' for multithreading)
| samtools sort \
    # Output in BAM format
    -O BAM \
    -T <temporary_file_directory> \
    - \
# Mark duplicate alignments from a coordinate sorted file that has been through 'samtools fixmate' with the '-m' option. This program relies on the MC and ms tags that fixmate provides. (can use '-@' for multithreading)
| samtools markdup \ 
    # Print some basic stats (optional)
    -s \
    # Name of stat file (optional)
    -f <markdup.stats> \
    -T <temporary_file_directory> \
    - \
    <complete_aligned_markdup.bam>
```

#### **Filtering**

```bash
# (can use '-@' for multithreading)
samtools view \
    # Output in BAM format
    -b \
    # Do not output alignments with any bits set in designated flag present in the FLAG field. '3844' excludes: unmapped reads, not primary alignment (alternative mappings when multiple mappings are presented), failed reads, optical or PCR duplicates, supplementary reads (the corresponding alignment is part of a chimeric alignment).
    -F 3844 \
    # Skip alignments with a MAPQ < 20
    -q 20 \
    -o <complete_aligned_filtered.bam> \
    <complete_aligned_markdup.bam>

# Index coordinated sorted BAM file for fast random access. (can use '-@' for multithreading)
samtools index \
    # Creates a BAI index. The BAI index format can handle individual chromosomes up to 512 Mbp in length.
    -b \
    <complete_aligned_filtered.bam> \
    > <complete_aligned_filtered.bam>.bai
```

#### **Extract unmapped reads**

```bash
# (can use '-@' for multithreading)
samtools view \
    # Output in BAM format
    -b \
    # Only output alignmnets that match with all bits set in designated flag present in the FLAG field. '4' includes: unmapped reads.
    -f 4 \
    -o <unmapped.bam> \
    <complete_aligned_markdup.bam>

# Index coordinated sorted BAM file for fast random access. (can use '-@' for multithreading)
samtools index \
    # Creates a BAI index. The BAI index format can handle individual chromosomes up to 512 Mbp in length.
    -b \
    <unmapped.bam> \
    > <unmapped.bam>.bai
```

#### **Post filtering stats**

```bash
# Change of JAVA_OPTS in qualimap script
shopt -s extglob
qualimap_path=$(dirname "$(which python)")/qualimap
grep -qxF '\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE"' "$qualimap_path" || sed -i '47s#.*#\tjava_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE"#' "$qualimap_path"

# Evaluate NGS mapping to a reference genome.
qualimap bamqc \
    # Input mapping file in BAM format
    -bam <complete_aligned_filtered.bam> \
    # Feature file with regions of interest in GFF/GTF format
    -gff <reference_genome.gff> | <reference_genome.gtf> \
    # Output folder for HTML report and raw data
    -outdir <qualimap/> \
    # Output file for PDF report
    -outfile <qualimap.pdf> \
    # Format of the output reprt. Default HTML.
    -outformat PDF \
    # Use this argument to set Java memory heap size.
    --java-mem-size=20G
```

### **03 - Initial Analysis Files, VCF and SFS**

-----

#### **Conda Enviroment**

Date: 01/06-2023

```yaml
name: vcf
channels:
  - genomedk
  - gwforg
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - freebayes=1.3.6
  - snpeff=5.1
  - snpsift=5.1
  - bioawk
  - ripgrep
  - bedtools=2.30.0
  - bcftools=1.16
  - vcftools=0.1.16
  - pandas
  - agat=1.0.0
  - rename=1.601
prefix: /home/jepe/miniconda3/envs/vcf
```

#### **Create VCF**

```bash
# Holds list of all SQ names from BAM file
SQlist=$(samtools view -H "$targetbam" | rg '@SQ' | awk '{print $2}')

# Turns SQlist into an array with SQ names sliced to easily work as individual region names
x=1
for line in $SQlist; do
    SQarray["$x"]="${line:3}"
    ((x++))
done

# Slices SQ names to work as region name. SQnum is the index number of the current region for analysis
region=${SQarray[SQnum]}

freebayes \
    # Use FILE as the reference sequence for analysis. An index file (FILE.fai) will be created if none exists. If neither --targets nor --region are specified, FreeBayes will analyze every position in this reference.
    -f <reference_genome> \
    # Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. (Set to 0 to use all). Default all.
    -n 3 \
    # Sets the default ploidy for the analysis to <ploidy>.
    -p <ploidy> \
    # Limit analysis to the specified region, 0-base coordinates, end_position not included (same as BED format). 'n' being array index for region
    -r "$region" \
    # Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position. Default '0.05'. If sample is individually resequenced use '0.05' if sample is pool resequenced use 0.
    --min-alternate-fraction 0 \
    # Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position.  Default '2'
    --min-alternate-count 2 \
    # Report even loci which appear to be monomorphic, and report all considered alleles, even those which are not in called genotypes. Loci which do not have any potential alternates have '.' for ALT.
    --report-monomorphic \
    # Assume that samples result from pooled sequencing. Model pooled samples using discrete genotypes across pools. When using this flag, set --ploidy to the number of alleles in each sample or use the --cnv-map to define per-sample ploidy.
    --pooled-discrete  \
    # File containing a list of BAM files to be analysed. Used if data is individually resequenced
    -L <bam_list> \
    # BAM file to be analysed. Used if data is pool resequenced
    <bam_file> \
| SnpSift intervals \
    # Exclude VCF entries in intervals
    -x \
    # BED file of repeat regions
    <repeat_regions> \
| bcftools \
    # Filter SNPs within 5 base pairs of an indel
    --SnpGap 5 \
| SnpSift filter \
    # Filters out any entry which doesn't have TYPE snp
    "(TYPE has 'snp')" \
    > "$region"_complete.vcf

# Creates empty filelist document and populated it with all regional VCF files
filelist=filelist.txt
echo -n "" > "$filelist"
for vcf in ./*_complete.vcf; do
    echo "$vcf" >> "$filelist"
done

# Concatenate VCF/BCF files. (can use '--threads' for multithreading)
bcftools concat \
    # Read the list of files from a file.
    -f "$filelist" \
    # Write output to a file
    -o <complete.vcf> \
    # Output type. v: uncompressed VCF
    -O v
```

#### **Create SnpEff Database**

```bash
# If annotation file is only available in GFF format, the following can convert it to GTF 2.2. GTF 2.2 is chosen as it is the newest GTF format compatible with SnpEff, and supposedly handles better than GFF format.
agat_convert_sp_gff2gtf.pl \
    # Input GFF/GTF file that will be read
    --gff <gff_file> \
    # Version of the GTF output. GTF2.2 (9 feature types accepted): CDS, start_codon, stop_codon, 5UTR, 3UTR, inter, inter_CNS, intron_CNS and exon. Default '3' (GTF3).
    --gtf_version 2.2 \
    # Output GTF file. If no output file is specified, the output will be written to STDOUT.
    -o <gtf_file>

# Assuming the use of a conda environment, the following will get path to your version of SnpEff and its data folder (if the data folder doesn't exist it is created)
snpeff_path=("$(dirname "$(dirname "$(which python)")")"/share/snpeff*)
snpeffdata="${snpeff_path[*]}"/data
[ -d "$snpeffdata" ] || mkdir "$snpeffdata"

# Creates backup of the original SnpEff configuration if there isn't one (Good to have a backup as changes will be made to the configuration file)
snpeffconfig_backup="${snpeff_path[*]}"/snpEff.config.backup
if [ ! -e "$snpeffconfig_backup" ]; then
    cp "${snpeff_path[*]}"/snpEff.config "$snpeffconfig_backup"
fi

# Add new entry to SnpEff configuration file
echo -e \
"\n# <species_name> genome, version <reference_genome>\n<reference_genome>.genome : <species_name>" \
>> "${snpeff_path[*]}"/snpEff.config
```

It should look something like this:

`# Orchesella cincta genome, version GCA_001718145.1_ASM171814v1`  
`GCA_001718145.1_ASM171814v1.genome : Orchesella cincta`

```bash
# Creates species folder within the SnpEff directory
speciesdata="$snpeffdata"/<reference_genome>
[ -d "$speciesdata" ] || mkdir "$speciesdata"

# Extract CDS sequences named using transcript ID and places it within the species folder
agat_sp_extract_sequences.pl \
    # Input GTF/GFF file.
    -g <gtf_file> \
    # Input fasta file.
    -f <reference_genome> \
    # Define the feature you want to extract the sequence from. Default 'cds'.
    -t cds \
    # Output fasta file. If no output file is specified, the output will be written to STDOUT.
    -o "$speciesdata"/cds.fa

# Extract protein sequences and places it within the species folder
agat_sp_extract_sequences.pl \
    # Input GTF/GFF file.
    -g <gtf_file> \
    # Input fasta file.
    -f <reference_genome> \
    # Define the feature you want to extract the sequence from. Default 'cds'.
    -t cds \
    # Will translate the extracted sequence in Amino acid.
    -p \
    # Output fasta file. If no output file is specified, the output will be written to STDOUT.
    -o "$speciesdata"/protein.fa

# Create copies of GTF file and genome reference sequence file in species folder renaming them to fit SnpEff naming convention
cp <gtf_file> "$speciesdata"/genes.gtf
cp <reference_genome> "$speciesdata"/sequences.fa

# Using explicit execution of snpEff.jar to change default java vm setting for heap size
java \
    # Sets initial heap size to 80G
    -Xms80G \
    # Sets maximum heap size to 80G
    -Xmx80G \
    # Build a SnpEff database.
    -jar "${snpeff_path[*]}"/snpEff.jar build \
        # Use GTF 2.2 format.
        -gtf22 \
        # Verbose mode. Makes potential troubleshooting easier
        -v \
        <reference_genome>
```

#### **Annotation with SnpEff**

```bash
# Assuming the use of a conda environment, the following will get path to your version of SnpEff and its configuration file
snpeff_path=("$(dirname "$(dirname "$(which python)")")"/share/snpeff*)
snpeffconfig="${snpeff_path[*]}"/snpEff.config

# Annotate variants / calculate effects
snpEff ann \
    # Create CSV summary file.
    -csvStats snpEff_summary.csv \
    # Specify config file.
    -c "$snpeffconfig" \
    # Verbose mode.
    -v \
    # The name given to the entry in the configuration file following following 'version'
    <genome_version> \
    # Input VCF file
    <complete.vcf> \
    > <complete.ann.vcf>
```

### **04 - Genome Assembly**

-----

#### **Conda Enviroment**

Date: 01/06-2023

Used for everything in the **BlobTools** section

```yaml
name: blobtools
channels:
  - genomedk
  - gwforg
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - blast=2.13.0
  - bioawk
  - minimap2=2.24
  - diamond=2.1.4
  - samtools=1.16.1
  - ripgrep
  - blobtools=1.1.1
prefix: /home/jepe/miniconda3/envs/blobtools
```

Date: 01/06-2023

Used for everything in the **BUSCO** section

```yaml
name: busco
channels:
  - genomedk
  - gwforg
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - busco=5.4.5
  - ripgrep
  - bioawk
prefix: /home/jepe/miniconda3/envs/busco
```

Date: 01/06-2023

Used for everything besides the **BlobTools** and **BUSCO** sections

```yaml
name: genome_assembly
channels:
  - genomedk
  - gwforg
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - 3d-dna=201008
  - bwa=0.7.17
  - pigz=2.6
  - bamtools=2.5.2
  - hifiasm=0.18.9
  - jellyfish=2.2.10
  - bioawk
  - matplotlib
  - purge_dups=1.2.6
  - blast=2.13.0
  - genomescope2=2.0
  - minimap2=2.24
  - pandas
prefix: /home/jepe/miniconda3/envs/genome_assembly
```

#### **HiFi adapter filtering**

```bash
# Convert .bam to .fastq and remove reads with remnant PacBio adapter sequences
# Adding script and database path to PATH variable
export PATH=$PATH:[PATH TO HiFiAdapterFilt]
export PATH=$PATH:[PATH TO HiFiAdapterFilt]/DB

# (can use '-t' for multithreading)
bash HiFiAdapterFilt/pbadapterfilt.sh \
    -p <PacBio_HiFi_sequence_file_prefix> \
    -o <output_directory_prefix>
```

<u>*Output files*</u>

- Output of BLAST search  
  `<PacBio_HiFi_sequence_file_prefix>.contaminant.blastout`  
- Headers of PacBio adapter contaminated reads to be removed  
  `<PacBio_HiFi_sequence_file_prefix>.blocklist`  
- FASTQ reads free of PacBio adapter sequence ready for assembly  
  `<PacBio_HiFi_sequence_file_prefix>.filt.fastq.gz`  
- File with simple math on number of reads reamoved etc.  
  `<PacBio_HiFi_sequence_file_prefix>.stats`  

#### **K-mer analysis**

```bash
# Count k-mers or qmers in FASTA or FASTQ files (can use '-t' for multithreading)
jellyfish count \
    # Count on both strands, canonical representation. Default 'false'. Needs to be set when running on sequencing reads. If running on an actual genome or finished sequence do NOT set.
    -C \
    -m <kmer_length> \
    # Hash size. Counting k-mers in short sequencing reads hash size should be approximately 'G + k * n)/0.8' where 'G' is the size of the genome, 'n' is the number of reads, and 'k' is the average number of errors per read. Counting k-mers in an assembled sequence the hash size should be 'G', which is the length of the assembled sequence. ISO suffixes can be used, thus '8M' stands for 8 million entries, while '8G' stands for 8 billion entries.
    -s 8G \
    # Default 'mer_counts'
    -o <output_prefix> \
    # If <sequence_file> is gzipped use process substitution with zcat.
    <sequence_file> | <(zcat <sequence_file>)
```

```bash
# Creates a histogram of k-mer occurrences (can use '-t' for multithreading)
jellyfish histo \
    # High count value of histogram. Default '10000'
    -h 100000 \
    -o <output_file> \
    <kmer_counts_file>
```

```bash
# Reference free profiling of polyploid genomes
genomescope2 \
    -i <kmer_histogram_file> \
    # To be sure, this will always be '2' in our case, as we're working with diploid animals
    -p <ploidy> \
    -o <output_directory> \
    # Should be the same as used by 'jellyfish count'
    -k <kmer_length> \
    -n <output_prefix> \
    # Optional maximum k-mer coverage threshold. K-mers with coverage greater than <max_kmer_coverage> are ignored by the model.
    -m <max_kmer_coverage>
```

#### **HiFi assembly**

```bash
# Fast haplotype-resolved de novo assembler initially designed for PacBio HiFi reads. (can use '-t' for multithreading)
hifiasm \
    # Prefix of output files. Default 'hifiasm.asm.'
    -o <output_prefix> \
    # Length that should be trimmed from ends of reads. Default '0'. Generally use default.
    -z <trim_length> \
    # Similarity threshold for duplicate haplotigs in read-level. Default '0.55' for '-l 3' If sample has high heterozygosity, then lower the value, start at 0.1.
    -s <similarity_threshold> \
    # Purge level. '0' for no purging, '1' for light purging, and '2'/'3' for aggressive purging. As we are dealing with very high levels of heterozygosity we usually always use '3'.
    -l 3 \
    # Output a primary assembly and an alternate assembly
    --primary \
    # Homozygous read coverage. Default 'auto'. Normally don't set this
    --hom-cov <INT> \
    <PacBio_HiFi_sequence_file>
```

<u>*Output files*</u>

- Assembly graph of primary contigs:  
  `<output_prefix>.bp.p_ctg.gfa`  
- Partially phased contig graph of haplotype1. Only gets produced when `--primary` is not set:  
  `<output_prefix>.bp.hap1.p_ctg.gfa`  
- Partially phased contig graph of haplotype2. Only gets produced when `--primary` is not set:  
  `<output_prefix>.bp.hap2.p_ctg.gfa`  
- Assembly graph of alternate contigs. Only gets produced when `--primary` is set:  
  `<output_prefix>.a_ctg.gfa`  
- Haplotype-resolved raw unitig graph. This graph keeps all haplotype information:  
  `<output_prefix>.r_utg.gfa`  
- Haplotype-resolved processed unitig graph without small bubbles. Small bubbles might be caused by somatic mutations or noise in data, which are not the real haplotype information. Hifiasm automatically pops such small bubbles based on coverage. The option `--hom-cov` affects the result. See homozygous coverage setting for more details. In addition, the option `-p` forcedly pops bubbles:  
  `<output_prefix>.p_utg.gfa`  

For each graph file, hifiasm also outputs a simplified version (`*noseq*gfa`) without sequences for the ease of visualization. The coordinates of low quality regions are written to `*lowQ.bed` in BED format.

#### **----- Making a FASTA formatted file from GFA**

```awk
awk '/^S/{print ">"$2"\n"$3}' <output_prefix>.bp.p_ctg.gfa | fold > <genome_assembly>.fasta
```

#### **BUSCO analysis**

```bash
# The Benchmarking Universal Single-Copy Ortholog assessment tool.
busco \
    # Force rewriting of existing files. Must be used when output files with the provided name already exists
    -f \
    -i <sequence_file> \
    # BUSCO analysis mode. 'genome' for genome assemblies (DNA). Can also be set to 'transcriptome' for transcriptome assemblies (DNA) or 'proteins' for annotated gene sets (protein)
    -m genome \
    # Output folders and files will be labelled with this name.
    -o <output_name> \
    # Location of results folder, excluding results folder name. Default is current working directory.
    --out_path <output_path> \
    # Local filepath for storing BUSCO dataset.
    --download_path /faststorage/project/EcoGenetics/BACKUP/database/busco/busco_downloads \
    # Specify the name of the BUSCO lineage to be used.
    -l /faststorage/project/EcoGenetics/BACKUP/database/busco/busco_downloads/lineages/arthropoda_odb10
```

#### **Purge dups**

When working with high heterozygosity assemblies several rounds of purge dup'ping should be utilized.

```bash
# Align PacBio data and generate PAF file. (can use '-t' for multithreading)
minimap2 \
    # Preset (always applied before other options). 'map-pb'/'map-ont' - PacBio CLR/Nanopore vs reference mapping, 'map-hifi' - PacBio HiFi reads vs reference mapping, 'ava-pb'/'ava-ont' - PacBio/Nanopore read overlap, 'asm5'/'asm10'/'asm20' - asm-to-ref mapping, for ~0.1/1/5% sequence divergence, 'splice'/'splice:hq' - long-read/Pacbio-CCS spliced alignment, 'sr' - genomic short-read mapping
    -x map-hifi \
    # File containing genome assembly corresponding to the selected raw PacBio HiFi reads.
    <genome_assembly> \
    # File containing raw PacBio HiFi reads.
    <PacBio_HiFi_sequence_file> \
# Gzip alignment to save space
| gzip \
    # Write on standard output, keep original files unchanged
    -c \
    - \
    # Alignment in gzipped PAF format
    > <alignment.paf.gz>
```

```bash
# Calculate read depth histogram and base-level read depth
pbcstat \
    -O <output_directory> \
    <alignment.paf.gz>
```

<u>*Output files*</u>

- `PB.stat`  
- `PB.base.cov`

```bash
# Calculate coverage cutoffs
calcuts \
    <PB.stat> \
    > <cutoffs>
```

```bash
# Split assembly
split_fa \
    <genome_assembly> \
    > <genome_assembly.split>
```

```bash
# Self-self alignment. (can use '-t' for multithreading)
minimap2 \
    # Assembly-to-reference mapping for ~0.1 sequence divergence
    -x asm5 \
    # If query sequence name/length are identical to the target name/length, ignore diagonal anchors. This option also reduces DP-based extension along the diagonal.
    -D \
    # Retain all chains and don’t attempt to set primary chains.
    -P \
    <genome_assembly.split> \
    <genome_assembly.split> \
# Gzip alignment to save space
| gzip \
    # Write on standard output, keep original files unchanged
    -c \
    - \
    # Alignment in gzipped PAF format
    > <self_alignment.paf.gz>
```

```bash
# Purge haplotigs and overlaps
purge_dups \
    # 2 rounds chaining. Default 'false'
    -2 \
    # Cutoffs file
    -T <cutoffs> \
    # Base level coverage file
    -c <PB.base.cov> \
    <self_alignment.paf.gz> \
    > <dups.bed>
```

```bash
# Get purged assembly and haplotig sequences from draft assembly
get_seqs \
    # This will only remove haplotypic duplications at the ends of the contigs. If you want to remove the duplications in the middle, remove the '-e' option at your own risk, as it may delete false positive. (We have run it without the '-e' option)
    -e \
    <dups.bed> \
    <genome_assembly>
```

<u>*Output files*</u>

- `purged.fa`

#### **BlobTools**

#### **----- COV file**

```bash
# (can use '-t' for multithreading)
minimap2 \
    # PacBio HiFi reads vs reference mapping
    -x map-hifi \
    # Output in SAM format
    -a \
    <genome_assembly> | <purged.fa> \
    <PacBio_HiFi_sequence_file> \
# Position sort alignment by leftmost coordinate. (can use '-@' for multithreading)
| samtools sort \
    -o <aligned.bam> \
    -O BAM \ Output in BAM format
    -T <temporary_file_directory> \
    -

# Index coordinated sorted BAM file for fast random access. (can use '-@' for multithreading)
samtools index \
    # Creates a BAI index. The BAI index format can handle individual chromosomes up to 512 Mbp in length.
    -b \
    <aligned.bam>
```

```bash
# COV is a custom file format written by BlobTools. It contains all the information needed for 'blobtools create'. Its usage is encourage as it is easier to parse.
blobtools map2cov \ 
    -i <genome_assembly> | <purged.fa> \
    -b <aligned.bam> \
    -o <output_prefix>
```

<u>*Output files*</u>

- `<output_prefix>.cov`

#### **----- HITS files**

Sequence Similarity search algorithms
Since the goal is to obtain hits for the sequences in the assembly, the queries of the sequence similarity searches are nucleotide sequences. The sequence similarity search algorithm determines which type (nucleotide or protein) of sequence collection can be used.  
The following sequence similarity search algorithms can be used:  

| Algorithm | Scope | Database | Speed |
|:-:|:-:|:-:|:-:|
| BLASTn megablast | Used to find very similar (e.g., intraspecies or closely related species) sequences | Nucleotide | Fastest |
| BLASTn dc-megablast | Used to find more distant (e.g., interspecies) sequences | Nucleotide | Slow |
| BLASTx | Used to find distant sequences | Protein | Slowest |
| Diamond blastx | Used to find distant sequences | Protein | Faster than BLASTx |

Recommendation from BlobTools documentation:  
Most accurate results are obtained with the following searches, supplied in this order and using taxrule ‘bestsumorder’.

The UniProt Reference Proteomes database is recommended when using DIAMOND blastx

```bash
# Align DNA query sequences against a protein reference database. (can use '--threads' for multithreading)
diamond blastx \
    # Directory to be used for temporary storage.
    --tempdir <temporary_file_directory> \
    # Path to the query input file in FASTA or FASTQ format (may be gzip compressed)
    --query <genome_assembly> | <purged.fa> \
    # The maximum number of target sequences per query to report alignments for. Default '25'
    --max-target-seqs 1 \
    # Enable sensitive mode designed for full sensitivity for hits of >40% identity.
    --sensitive \
    # Path to DIAMOND database file
    --db /faststorage/project/EcoGenetics/BACKUP/database/uniprot/uniprot_ref_proteomes.dmnd \
    # Maximum e-value to report alignments. Default '0.001'
    --evalue 1e-25 \
    # Format of the output file. '0' is BLAST pairwise, '5' is BLAST XML, '6' is BLAST tabular, '100' is DIAMOND alignment archive (DAA), '101' is SAM
    --outfmt 6 \
    # Path to output file
    --out <diamondblast.out>

# Using the BlobTools module taxify, the user can create their own Hits file. It must be: based on other TSV input, such as the output of other contaminant screening pipelines, and based on a single NCBI TaxID and a score
blobtools taxify \
    # BLAST/DIAMOND similarity search result (TSV format).
    -f <diamondblast.out> \
    # TaxID mapping file (contains seqid and taxid)
    -m /faststorage/project/EcoGenetics/BACKUP/database/uniprot/uniprot_ref_proteomes.taxids \
    # Zero-based column of sseqid in TaxID mapping file (it will search for sseqid in this column)
    -s 0 \
    # Zero-based Column of taxid in TaxID mapping file (it will extract for taxid from this column)
    -t 2 \
    -o <output_directory_prefix>
```

<u>*Output files*</u>

- `<diamondblast>.taxified.out`

The NCBI nt database is recommend when using BLASTn

```bash
# (can use '-num_threads' for multithreading)
blastn \
    # Traditional megablast used to find very similar (e.g., intraspecies or closely related species) sequences
    -task megablast \
    # Input file
    -query <genome_assembly> | <purged.fa> \
    # Path to BLAST database
    -db /faststorage/project/EcoGenetics/BACKUP/database/blastdb/nt \
    # Outfmt: 6 = tabular, additional formatting, std = standard. Resulting format: qseqid (Query Seq-ID), staxids (Unique Subject Taxonomy IDs, separated by a ';'), bitscore (Bit Score), pident (Percentage of identical matches), length (Alignment lenght), mismatch (Number of mismatches), gapopen (Number of gap openings), qstart (Start of alignment in query), qend (End of alignment in query), sstart (Start of alignment in subject), send (End of alignment in subject), evalue (Expect value)
    -outfmt '6 qseqid staxids bitscore std' \
    # Maximum number of aligned sequences to keep. Default '500'.
    -max_target_seqs 10 \
    # Set maximum number of HSPs per subject sequence to save for each query
    -max_hsps 1 \
    # Expectation value (E) threshold for saving hits. Default '10'
    -evalue 1e-25 \
    -out <ncbimegablast.out>

```

#### **----- Blobplot**

```bash
blobtools create \
    # NCBI nodes.dmp file.
    --nodes /faststorage/project/EcoGenetics/BACKUP/database/ncbi_taxonomy/nodes.dmp \
    # NCBI names.dmp file.
    --names /faststorage/project/EcoGenetics/BACKUP/database/ncbi_taxonomy/names.dmp \
    # FASTA file of assembly.
    -i <genome_assembly> | <purged.fa> \
    # Hits file in format (qseqid\ttaxid\tbitscore) (e.g. BLAST output "--outfmt '6 qseqid staxids bitscore'"). Can be specified multiple times
    -t <diamondblast>.taxified.out \
    -t <ncbimegablast.out> \
    # Taxrule determines how taxonomy of blobs is computed (by default both are calculated). 'bestsum': sum bitscore across all hits for each taxonomic rank, 'bestsumorder': sum bitscore across all hits for each taxonomic rank (If first <TAX> file supplies hits, bestsum is calculated. If no hit is found, the next <TAX> file is used.)
    -x bestsumorder \
    # COV file(s), can be specified multiple times
    -c <COV> \
    # BlobDB output prefix
    -o <output_prefix>
```

```bash
blobtools plot \
    # BlobDB file (created with "blobtools create")
    -i <blobdb>.blobDB.json \
    -o <output_directory_prefix> \
    # Taxrule which has been used for computing taxonomy (Supported: 'bestsum', 'bestsumorder'). Default 'bestsum'
    -x bestsumorder
```

<u>*Output files*</u>

- BlobPlots are two-dimensional scatter plots, decorated with coverage and GC histograms. Sequences are represented by circles (with diameters proportional to sequence length) in the scatter plot and coloured by taxonomic affiliation. Circles are positioned on the Y-axis in the scatter plot based on the base coverage of the sequence in the coverage library, a proxy for molarity of input DNA in the sequencing reaction. Circles are position on the X-axis based on their GC content, the proportion of G and C bases in the sequence, which can differ substantially between genomes. Coverage and GC histograms are drawn for each taxonomic group, which are weighted by the total span (cumulative length) of sequences occupying each bin. A legend reflects the taxonomic affiliation of sequences and lists count, total span and N50 by taxonomic group. Taxonomic groups can be plotted at any taxonomic rank and colours are selected dynamically from a colour map. The number of taxonomic groups to be plotted can be controlled (--plotgroups, default is ‘7’) and remaining groups are binned into the category ‘others’  
  `<blobdb>.blobDB.jason.bestsumorder.phylum.p8.span.100.blobplot.cov0.png`  
- ReadCovPlots visualise the proportion of reads of a library that are unmapped or mapped, showing the percentage of mapped reads by taxonomic group, as barcharts. ReadCovPlots can be of use for rapid taxonomic screening of multiple sequencing libraries within a single project. The underlying data of ReadCovPlots and additional metrics are written to tabular text files for custom analyses by the user.  
  `<blobdb>.blobDB.jason.bestsumorder.phylum.p8.span.100.blobplot.read_cov0.png`
- `<blobdb>.blobDB.jason.bestsumorder.phylum.p8.span.100.blobplot.stats.txt`

```bash
# Using 'blobtools view', information stored in a blobDB can be extracted. Using the default settings, blobtools view will generate a tabular output for the taxonomic rank of "phylum".
blobtools view \
    # BlobDB file (created with "blobtools create")
    -i <blobdb>.blobDB.json \
    -o <output_directory_prefix> \
    # Taxrule which has been used for computing taxonomy (Supported: 'bestsum', 'bestsumorder'). Default 'bestsum'
    -x bestsumorder \
    # Taxonomic rank(s) at which output will be written. (supported: 'species', 'genus', 'family', 'order',phylum', 'superkingdom', 'all'). Default 'phylum'
    -r all
```

<u>*Output files*</u>

| Column header | Description |
|:-:|:-:|
| Name | Name of the sequence |
| Length | Total length of the sequence, i.e. count(A, G, C, T, N) |
| GC | GC content percentage of the sequence, i.e. count(G, C)/count(A, G, C, T) |
| N | Number of N's in the sequence, i.e. count(N) |

- `<blobdb>.blobDB.table.txt`

### **05 - Repeat Detection, Annotating and Masking**

-----

Date: 01/06-2023

```yaml
name: repeatmasking
channels:
  - genomedk
  - gwforg
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - ripgrep
  - seqkit=2.3.1
  - repeatmasker=4.1.2.p1
  - bedtools=2.30.0
  - repeatmodeler=2.0.3
  - bioawk
prefix: /home/jepe/miniconda3/envs/repeatmasking
```

#### **Building database for RepeatModeler2**

```bash
# Format FASTA files for use with RepeatModeler
BuildDatabase \
    # The name of the database to create. We would typically use the species name
    -name <db_name> \
    # The name of the search engine we are using. Supported: 'abblast'/'wublast' or 'rmblast'.
    -engine rmblast \
    <genome_assembly>
```

#### **Creating model with RepeatModeler2**

Creates semi-standard species specific ID code. Ex. turns *Aglais_urticae* into *AglUrt*, followed by the number *1*, indicating that it is version 1. It is assumed that species name is written with '_' instead of the 'space' character.

```bash
genus="${<species_name>%_*}"
genus="${genus::3}"
genus="${genus^}"
species="${species_name#*_}"
species="${species::3}"
species="${species^}"
species_code="$genus""$species"1

# Model repetitive DNA
RepeatModeler \
    # Path to database created with 'BuildDatabase'
    -database <db_name> \
    # Run the LTR structural discovery pipeline and combine results with the RepeatScout/RECON pipeline
    -LTRStruct \
    # Specify the number of parallel search jobs to run. RMBlast jobs will use 4 cores each and ABBlast jobs will use a single core each. i.e. on a machine with 12 cores and running with RMBlast you would use -pa 3 to fully utilize the machine. (I have arbitrarily chosen '8' as it equals 32 cores)
    -pa 8
```

<u>*Output files*</u>

- All RepeatModeler2 output is normally placed into a directory name `RM_<datainfo>`. Important output is the final, classified repeat consensus libraries in FASTA format:  
  `<db_name>-families.fa`

Adds species code to outputted fasta file

```bash
cat <db_name>-families.fa \
| seqkit fx2tab \
| awk -v species_code=<species_code> '{ print species_code"_"$0 }' \
| seqkit tab2fx \
    > <db_name>-families.prefix.fa
```

Split fasta file into classified and unclassified elements

```bash
cat <db_name>-families.prefix.fa \
| seqkit fx2tab \
| rg -v "Unknown" \
| seqkit tab2fx \
    > <db_name>-families.prefix.fa.known
cat <db_name>-families.prefix.fa \
| seqkit fx2tab \
| rg "Unknown" \
| seqkit tab2fx \
    > <db_name>-families.prefix.fa.unknown
```

#### **Repeat masking based on GIRI RepBase**

!MAKE SETUP FOR DOWNLOADING AND SETTING UP GIRI REPBASE FOR ARTHROPODA!

```bash
# Mask repetitive DNA
RepeatMasker \
    # Use an alternate search engine to the default. Note: 'ncbi' and 'rmblast' are both aliases for the rmblastn search engine engine. The generic NCBI blastn program is not sensitive enough for use with RepeatMasker at this time. Supported: 'crossmatch', 'wublast', 'abblast', 'ncbi', 'rmblast', 'hmmer'
    -e rmblast \
    # The number of sequence batch jobs [50kb minimum] to run in parallel. RepeatMasker will fork off this number of parallel jobs, each running the search engine specified. For each search engine invocation ( where applicable ) a fixed the number of cores/threads is used: 'RMBlast': 4 cores, 'ABBlast': 4 cores, 'nhmmer': 2 cores, 'crossmatch': 1 core. To estimate the number of cores a RepeatMasker run will use simply multiply the -pa value by the number of cores the particular search engine will use.
    -pa 8 \ 
    # Writes output to this directory (default is query file directory, "-dir ." will write to current directory).
    -dir <output_directory> \
    # Returns repetitive regions in lowercase (rest capitals) rather than masked (This is called 'soft masking')
    -xsmall \
    # Allows use of a custom library (e.g. from another species)
    -lib /faststorage/project/EcoGenetics/BACKUP/database/giri_repbase/RepBaseCustom15.02.fasta/Arthropoda.rep \
    <genome_assembly>
```

<u>*Output files*</u>

- Full RepeatMasker results (file can get large)  
  `<genome_assembly>.cat.gz`  
- The masked reference genome. So, a copy of the original reference genomes with repeat elements soft masked  
  `<genome_assembly>.masked`  
- A tabular list of all repeat regions/elements annotated in this round of annotation/masking  
  `<genome_assembly>.out`  
- A summary table of the overall composition of the genome based on major repeat element groupings  
  `<genome_assembly>.tbl`  

Though `RepeatMasker` has the option to output a `GFF` formatted file of the repeat regions it is in an outdated `GFF2` format. The custom `AWK` script below uses the main output file from RepeatMasker to create a `GFF3` formatted file of the repeat regions:

```awk
# Make GFF 3 of repeat regions from RepeatMasker output file
awk \
    -F " " \
    'BEGIN {OFS="\t"; print "##gff-version 3"}
    NR > 3
        {if ($9 == "C")
            {strand = "-"}
        else
            {strand = "+"}
        if ($12 ~ /\(/)
            {start = $14}
        else {start = $12}
        print $5, "RepeatMasker", "repeat_region", $6, $7, ".", strand, ".", "ID="$15";Name="$10";Class="$11";Family="$11";Target="$10" "start" "$13}' \
    <genome_assembly>.out \
    > <genome_assembly>.repeats_giri.gff
```

#### **Repeat masking based on RepeatModeler2 model**

This round is run ontop of the resulting soft masked sequence from the previous round

```bash
RepeatMasker \
    # Use an alternate search engine to the default. Note: 'ncbi' and 'rmblast' are both aliases for the rmblastn search engine engine. The generic NCBI blastn program is not sensitive enough for use with RepeatMasker at this time. Supported: 'crossmatch', 'wublast', 'abblast', 'ncbi', 'rmblast', 'hmmer'
    -e rmblast \
    # The number of sequence batch jobs [50kb minimum] to run in parallel. RepeatMasker will fork off this number of parallel jobs, each running the search engine specified. For each search engine invocation ( where applicable ) a fixed the number of cores/threads is used: 'RMBlast': 4 cores, 'ABBlast': 4 cores, 'nhmmer': 2 cores, 'crossmatch': 1 core. To estimate the number of cores a RepeatMasker run will use simply multiply the -pa value by the number of cores the particular search engine will use.
    -pa 8 \
    # Writes output to this directory (default is query file directory, "-dir ." will write to current directory).
    -dir <output_directory> \
    # Returns repetitive regions in lowercase (rest capitals) rather than masked (This is called 'soft masking')
    -xsmall \
    # Allows use of a custom library (e.g. from another species)
    -lib <db_name>-families.prefix.fa \
    <genome_assembly>.masked
```

Small renaming of output files to avoid the naming convention 'masked.masked'. Otherwise output files are the same as in the last round.

```bash
rename 's/.masked././' ./*
```

Though `RepeatMasker` has the option to output a `GFF` formatted file of the repeat regions it is in an outdated `GFF2` format. The custom `AWK` script below uses the main output file from RepeatMasker to create a `GFF3` formatted file of the repeat regions:

```awk
awk \
    -F " " \
    'BEGIN {OFS="\t"; print "##gff-version 3"}
    NR > 3
        {if ($9 == "C")
            {strand = "-"}
        else
            {strand = "+"}
        if ($12 ~ /\(/)
            {start = $14}
        else {start = $12}
        print $5, "RepeatMasker", "repeat_region", $6, $7, ".", strand, ".", "ID="$15";Name="$10";Class="$11";Family="$11";Target="$10" "start" "$13}' \
    <genome_assembly>.out \
    > <genome_assembly>.repeats_rm.gff
```

#### **Combining results from repeat masking**

Combine full RepeatMasker result files

```bash
cat \
    # From the run on GIRI RepBase
    <genome_assembly>.cat.gz \
    # From the run on the RepeatModeler2 model
    <genome_assembly>.cat.gz \
    > <genome_assembly>.full_mask.cat.gz
```

Combine RepeatMasker tabular files

```bash
cat \
    # From the run on GIRI RepBase
    <genome_assembly>.out \
    <(tail \
        # Output the last 'K' lines, instead of the last 10; or use '-n +K' to output starting with the 'K'th line. The header, which is identical in the two files, covers the first 3 lines.
        -n +4 \
        # From the run on the RepeatModeler2 model
        <genome_assembly>.out \
    > <genome_assembly>.full_mask.out
```

To make combined versions of `<genome_assembly>.tbl` we use `ProcessRepeats` to resummarize repeat compositions from combined analysis of both RepeatMasker rounds

```bash
# Post process results from RepeatMasker and produce an annotation file.
ProcessRepeats \
    # Allows use of a custom library (e.g. from another species)
    -lib /faststorage/project/EcoGenetics/BACKUP/database/giri_repbase/RepBaseCustom15.02.fasta/Arthropoda.rep \
    <genome_assembly>.full_mask.cat.gz
```

Though `ProcessRepeats` has the option to output a `GFF` formatted file of the repeat regions it is in an outdated `GFF2` format. The custom `AWK` script below uses the main output file from RepeatMasker to create a `GFF3` formatted file of the repeat regions:

```awk
awk \
    -F " " \
    'BEGIN {OFS="\t"; print "##gff-version 3"}
    NR > 3
        {if ($9 == "C")
            {strand = "-"}
        else
            {strand = "+"}
        if ($12 ~ /\(/)
            {start = $14}
        else {start = $12}
        print $5, "RepeatMasker", "repeat_region", $6, $7, ".", strand, ".", "ID="$15";Name="$10";Class="$11";Family="$11";Target="$10" "start" "$13}' \
    <genome_assembly>.full_mask.out \
    > <genome_assembly>.repeats.gff
```

#### **Complete soft masking of genome assembly and create BED file of repeat regions**

```bash
# Mask a fasta file based on feature coordinates.
bedtools maskfasta \
    # Enforce "soft" masking. Mask with lower-case bases, instead of masking with Ns.
    -soft \
    # Input FASTA file
    -fi <genome_assembly> \
    # BED/GFF/VCF file of ranges to mask in '-fi'
    -bed <genome_assembly>.repeats.gff \
    # Output FASTA file
    -fo <species_code>.full_mask.soft.fna
```

Make *`BED3`* file of repeat regions. As GFF files report position in 1-based coordinates 1 is subtracted from the starting position to make it 0-based, which is the standard for `BED` format

```awk
awk \
    -F "\t" \
    'BEGIN {OFS = "\t"}
    {if ($0 ~ /^[^#]/)
        {print $1, ($4 - 1), $5}
    }' \
    <genome_assembly>.repeats.gff \
    > <species_code>.repeats.bed
```

### **Miscellaneous**

-----

#### **Make DIAMOND database from UniProt Reference Proteomes**

```bash
# Took 5 days
wget -q ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Reference_Proteomes_2022_05.tar.gz

# Unpack archive
for file in ./*.tar.gz; do
    tar -zxf "$file" -C .
done

# Unpack protein FASTAs for each kingdom
for file in ./*/*/*.fasta.gz; do
    if [[ ! "$file" == *DNA.fasta.gz ]];then
        if [[ ! "$file" == *additional.fasta.gz ]]; then
            gzip -d "$file"
        fi
    fi
done

for file in ./*/*/*.idmapping.gz; do
    gzip -d "$file"
done

# Concatenate all protein sequences into 'uniprot_ref_proteomes.fasta'
cat ./*/*/*.fasta > uniprot_ref_proteomes.fasta

# Simplyfy sequence IDs
cat uniprot_ref_proteomes.fasta | sed -r 's/(^>sp\|)|(^>tr\|)/>/g' | cut -f1 -d"|" > temp; mv temp uniprot_ref_proteomes.fasta

# Subset mapping files to only contain NCBI TaxID entries
cat ./*/*/*.idmapping | rg -j 2 "NCBI_TaxID" > uniprot_ref_proteomes.taxids

# Make Diamond DB
diamond makedb \
    --threads 2 \
    --in uniprot_ref_proteomes.fasta \
    -d uniprot_ref_proteomes
```

-----

#### **Setting up Progressive cactus aligner**

First create and activate a conda environment with 'virtualenv' installed.
Then cactus needs to be downloaded and unpacked in a desired directory:

```bash
# Download most recent version of cactus pre-compiled binaries
wget -q https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.6.0/cactus-bin-v2.6.0.tar.gz
# Decompress and unpack tarball
tar -xzf cactus-bin-v2.6.0.tar.gz
# Change working directory
cd cactus-bin-v2.6.0
```

Cactus will be run through its own virtual environment. This is setup and prepared the following way:

```bash
# Create virtual environment
virtualenv -p python 3.9 cactus_env
# Add PATH variable changes to 'activate' script
echo "export PATH=$(pwd)/bin:\$PATH" >> cactus_env/bin/activate
echo "export PYTHONPATH=$(pwd)/lib:\$PYTHONPATH" >> cactus_env/bin/activate
# Activate virtual envrionment
source cactus_env/bin/activate
# Install/upgrade series of dependencies listed within cactus directory using pip
python3 -m pip install -U setuptools pip
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt
```

Some tools are needed that are not included in the cactus directory. These are downloaded and execution permission set with the following:

```bash
cd bin
for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed axtChain pslPosTarget bedSort hgGcPercent mafToBigMaf hgLoadMafSummary; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/${i}; chmod +x ${i}; done
```
