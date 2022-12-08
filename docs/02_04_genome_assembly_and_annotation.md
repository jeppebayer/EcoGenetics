# GENOME ASSEMBLY AND ANNOTATION

[Previous](02_03_initial_analysis_procedure.md) | [Next](03_00_terminology.md)

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

## Genome Assembly

![jilong_assembly](../resources/jilong_assembly.png)

Genome assembly, including juicer aligning Hi-C data, scaffolding 3D-DNA, and exporting the scaffolded genome fasta file after manual review.

[Assembly_worklow](../scripts/jilong/genome_assembly_workflow.py)

Tips/tricks from Jilong to keep in mind:

- Remember to "wrap" the fasta file to be scaffolded when exporting fasta, it makes the process faster. By wrap, it means the fasta file is truncated into a fixed maximum number of characters per line. 3D-DNA come together with an embedded script to do this
- For me, running 3D-DNA without mis-joint correction rounds works well. (By setting -r 0). The mis-joint corrections sometimes may even decrease the assembly quality. It is because the combination of HiFi contigs and HiC provide quite high confidence already and mis-joint correction focuses more on troubles from short reads contigs.
- It is the XXXX.0.assembly file that i used to check HiC, fix orders, and split chromosomes after scaffolding.

## Genome Annotation

![jilong_annotation1](../resources/jilong_annotation1.png)
![kilong_annotation2](../resources/jilong_annotation2.png)

BRAKER2 has conda distribution and it is way easier to install it through conda than to manually install it

[Annotation_workflow](../scripts/jilong/genome_annotation_workflow.py)

Tips/tricks from Jilong to keep in mind:

- The repbase (curated repeat sequence database) is stored in the cluster /home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/repbase
- I use the repbase combined with repeatmodeler generated database to do repeatmasking. But in my case, I didn't see a huge difference in using repbase or not

[Previous](02_03_initial_analysis_procedure.md) | [Next](03_00_terminology.md)