#!/bin/env python
import sys

def parse_fasta(fasta_file: str):
    """Parses `FASTA` file returning all sequence names and lengths paired in a list of dictionaries.
    
    ::
    
        return [{'sequence_name': str, 'sequence_length': int}, ...]
    
    :param str fasta_file:
        Sequence file in `FASTA` format.
    """
    fasta_list = []
    seq_name = None
    length = 0
    with open(fasta_file, 'r') as fasta:
        for entry in fasta:
            entry = entry.strip()
            if entry.startswith(">"):
                if seq_name:
                    fasta_list.append({'sequence_name': seq_name, 'sequence_length': length})
                    length = 0
                entry = entry.split(" ", 1)
                seq_name = entry[0][1:]
            else:
                length += len(entry)
        fasta_list.append({'sequence_name': seq_name, 'sequence_length': length})
    return fasta_list

# Gets fasta_file parameter from argument 1. If second argument is given uses this as output file name, otherwise writes to stdout
fasta_file = sys.argv[1]
if sys.argv[2]:
    chrom_sizes = sys.argv[2]
else:
    chrom_sizes = sys.stdout

# Runs parse_fasta to parse given fasta file and output each entry as a string with tab-separated sequence name and length
with open(chrom_sizes, 'w') as sizes:
    for entry in parse_fasta(fasta_file):
        sizes.write('{name}\t{length}\n'.format(name=entry['sequence_name'], length=entry['sequence_length']))