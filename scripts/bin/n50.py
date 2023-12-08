#!/bin/env python
import os, sys, gzip

usage = f"\nUsage: {sys.argv[0]} [input file] [output file]\n"

def load_data(x: str):
    """Loads data either from standard input or from argument position 1. Input file can be gzipped.
    
    :param any x:
        Input data or '-'."""
    if x == '-':
        data = sys.stdin
    elif x.endswith(".gz"):
        data = gzip.open(x, 'rt')
    else:
        data = open(x, 'r')
    return data

def n50(lengths: str | int):
    """Calculates total sequence length and N50 for a series of sequence lengths.
    
    :param str lengths:
        Series on sequence lengths."""
    cumulative_sum = 0
    n = 0
    rows = (row for row in lengths)
    sequence_lengths = sorted([int(entry.rstrip().split("\t")[1]) for entry in rows], reverse = True)  
    total_sequence_length = sum(sequence_lengths)
    for sequence_length in sequence_lengths:
        cumulative_sum += sequence_length
        n += 1
        if cumulative_sum >= (total_sequence_length // 2):
            yield f'Total length\t{total_sequence_length}\nn\t{len(sequence_lengths)}\nN50\t{sequence_length}\nn\t{n}'
            break

def continuously_write_output(input: str, output: str):
    """Continuously writes each iteration of input to output"""
    with output as outfile:
        for entry in input:
            outfile.write(f"{entry}\n")
            outfile.flush()

if len(sys.argv) <= 1 or len(sys.argv) >= 4:
    sys.stdout.write(f'{usage}\n')
    exit(1)
elif len(sys.argv) == 2:
    continuously_write_output(n50(load_data(sys.argv[1])), sys.stdout)
else:
    if os.path.exists(sys.argv[2]):
        os.remove(sys.argv[2])
    continuously_write_output(n50(load_data(sys.argv[1])), open(sys.argv[2], 'w'))