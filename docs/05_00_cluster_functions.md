# **CLUSTER FUNCTIONS**

[Previous](04_00_software.md)

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

## **Request interactive computing node:**

```json
srun --mem-per-cpu=16g --time=3:00:00 --account=EcoGenetics --pty bash
```

## **Request Jupyter node**

From own machine in popgen environment

```json
slurm-jupyter -u jepe -A EcoGenetics -e jupyter -m 8g -t 3h --run notebook
```

## **Change permission for files and folders:**

```bash
chmod [options] [reference][operator][modes] file
```

[reference]: `u` (user), `g` (group), `o` (others), `a` (all)  
[operator]: `+` (add permissions), `-` (remove permissions)  
[mode]: `r` (read), `w` (write), `x` (execute)  
[options]: `-R` (recursive, i.e. include all objects in subdirectories)

Example:

```bash
chmod g+rwx [filename] 
```

Add read, write and execute permissions to all group users

## **Interact with the queue:**

Interaction with the queue is done using BASH script  
First lines of BASH file:

```bash
#!/bin/bash
#SBATCH [flag1]
#SBATCH [flag2]
```

Common flags:

- `--account` Account to submit the job under.
- `--partition`         One or more comma-separated partitions that the job may run on. Jobs submitted to the gpu partition should also use the –gres flag.
- `--mem-per-cpu`       Memory allocated per allocated CPU core.
- `--cpus-per-task`     Number of cores allocated for the job. All cores will be on the same node.
- `--ntasks`            Number of cores allocated for the job. Cores may be allocated on different nodes.
- `--nodes`             Number of nodes allocated for the job. Can be combined with -n and -c.
- `--time`              Maximum time the job will be allowed to run.
- `--constraint`        Constrain nodes to be allocated.

Send to queue:

```bash
sbatch [filename]
```

Get status on job:

```bash
jobinfo [jobnumber]
```

Cancel job:

```bash
scancel [jobnumber]
```

Current jobs in queue by user:

```bash
squeue -u [usersname]
```

## Send files to/from the cluster

```bash
scp ./file jepe@login.genome.au.dk:dir/
```

```bash
scp jepe@login.genome.au.dk:dir/file .
```

[Previous](04_00_software.md)