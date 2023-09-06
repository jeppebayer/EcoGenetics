# Progressive Cactus

Cactus takes a set of **softmasked** genome assemblies in `FASTA` format (can be gzipped) and a phylogenetic tree realting them as input.
Each genome should correspond to a leaf in the tree. Cactus outputs a multiple genome alignment in `HAL` (Hierachical Alignment) format, which includes ancestral sequences, one for each internal node in the tree.

The basic format to run Cactus is:

```bash
cactus <JobStorePath> <SeqFile> <OutputHAL>
```

- `JobStorePath` is where intermediate files, as well as job metadata is kept by Toil (Toil is an open-source pure-Python workflow engine which is used by Cactus). This path must be accessible to all worker systems.
- `SeqFile` is a text file containing the locations of the input sequences as well as their phylogenetic tree. The tree will be used to progressively decompose the alignment by iteratively aligning sibling genomes to estimate their parents in a bottom-up fashion. Polytomies in the tree are allowed, though the amount of computation required for a sub-alignment rises quadratically with the degree of the polytomy. Cactus uses the predicted branch lengths from the tree to determine appropriate pairwise alignment parameters, allowing closely related species to be aligned more quickly with no loss in accuracy.  
An optinal '\*' can be placed at the beginning of a name to specify that its assembly is of reference quality. This implies that it ca nbe used as an outgroup for sub-alignments. If no genomes are marked in this way, all genomes are assumed to be of reference quality. The '\*' should only be placed on the name-path lines nad not inside the tree.  
The `SeqFile` must conform to the following guidelines:

  - The tree must be on a single line. All leaves must be labeled and these labels must be unique. Ancestors may be named, or left blank (in which case the ancestors in the final output will automatically be labeled Anc0, Anc1, etc.) Labels must not contain any spaces.
  - Branch lengths that are not specified are assumed to be 1.
  - Lines beginning with # are ignored.
  - Sequence paths must point to either a `FASTA` file or a directory containing 1 or more `FASTA` files.
  - `FASTA` files may be gzipped.
  - Sequence paths must not contain spaces.
  - Each name / path pair must be on its own line.
  - `http://`, `s3://`, etc. URLs may be used.

  The format of the `SeqFile` must be as follows:  
  
        NEWICK tree
        name1 path1
        name2 path2
        ...
        nameN pathN

    **Ensure that genomes are soft-masked** (ideally with RepeatMasker). Genomes that are not properly masked can take tens of times longer to align than those that are masked. **Hard-masking is strongly discouraged**.
- `OutputHAL` represents the multiple alignment, including all input and inferred ancestral sequences. It is stored in `HAL` format, and can be accessed with HAL tools, which are all included in Cactus.

<https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md>