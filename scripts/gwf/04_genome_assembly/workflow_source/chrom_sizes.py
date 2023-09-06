#!/bin/env python
import sys

fasta_file = sys.argv[1]
if sys.argv[2]:
    chrom_sizes = sys.argv[2]
else:
    chrom_sizes = sys.stdout

def parse_fasta(fasta_file):
    """Function to parse FASTA file, retrieving all sequence names and lengths and returning them paired in a list of dictionaries."""
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

with open(chrom_sizes, 'w') as sizes:
    for entry in parse_fasta(fasta_file):
        sizes.write('{name}\t{length}\n'.format(name=entry['sequence_name'], length=entry['sequence_length']))