# Cluster commands

## Interactive job

You can request an interactive node on the cluster, which can be useful if you just testing you scripts

```bash
srun --cpus-per-task=<num_cpus> --mem=<ram>g --time=<duration> --account=EcoGenetics --pty bash
```

Here you can specify the number of cpus you want in \<num_cpus>, memory in GB ram in \<ram> and how long you would like access with \<duration> which is written in the format HH:MM:SS.

## Change permissions

```bash
chmod [options] [reference][operator][mode] [file | directory]
```

- \[reference\]: `u` (user), `g` (group), `o` (others), `a` (all)  
- \[operator\]: `+` (add permissions), `-` (remove permissions)  
- \[mode\]: `r` (read), `w` (write), `x` (execute)
- \[options\]: `-R` (recursive, i.e. include all objects in subdirectories)

Example:

```bash
chmod g+rwx file.txt
```

## Interacting with the queue

Interaction with the queue on the cluster is done using `bash` scripts.
The first lines of any `bash` script sent to the queue should be:

```bash
#!/bin/bash
#SBATCH [flag1]
#SBATCH [flag2]
...
#SBATCH [flagn]
```

Common flags:

- `--account`: Account to submit the job under.
- `--mem`: Memory allocated per allocated CPU core.
- `--cpus-per-task`: Number of cores allocated for the job. All cores will be on the same node.
- `--time`: Maximum time the job will be allowed to run.

For instance if I'm part of the *EcoGenetics* project and I want my job to run on *4 cores*, have access to *40 GB of RAM* shared on the cores and have reserve *6 hours* for my job to run, I would start my `bash` script with:

```bash
#!/bin/bash
#SBATCH --account=EcoGenetics
#SBATCH --cpus-per-task=4
#SBATCH --mem=40g
#SBATCH --time=06:00:00
```

