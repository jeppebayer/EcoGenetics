#!/bin/env python
import os, sys, gzip

usage = f"\nUsage: {sys.argv[0]} fa|fq [input file] [output file]\n"

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

def seq_length(sequence_file: str, filetype: str):
    """Get length of each sequence in :format:`fasta` and :format:`fastq` files.
    
    :param str sequence_file:
        Input file.
    :param str filetype:
        Filetype, either 'fa' (:format:`fasta`) or 'fq' (:format:`fastq`)."""
    if filetype not in ['fa', 'fq']:
        print(f'Unrecognized filetype: \'{filetype}\'')
        exit(1)
    rows = (row for row in sequence_file)
    entries = (entry.rstrip() for entry in rows)
    if filetype == 'fa':
        seq_name = None
        length = 0
        for entry in entries:
            if entry.startswith('>'):
                if seq_name:
                    yield f"{seq_name}\t{length}"
                    length = 0
                seq_name = entry.split(" ", 1)[0][1:]
            else:
                length += len(entry)
        yield f"{seq_name}\t{length}"
    if filetype == 'fq':
        for entry in entries:
            if entry.startswith('@'):
                yield f"{entry[1:]}\t{len(next(entries))}"
                next(entries)
                next(entries)

def continuously_write_output(input: str, output: str):
    """Continuously writes each iteration of input to output"""
    with output as outfile:
        for entry in input:
            outfile.write(f"{entry}\n")
            outfile.flush()

if len(sys.argv) <= 2 or len(sys.argv) >= 5:
    sys.stdout.write('{}\n'.format(usage))
    exit(1)
elif len(sys.argv) == 3:
    if sys.argv[1].endswith('.fa') or sys.argv[1].endswith('.fasta'):
        filteype = 'fa'
    elif sys.argv[1].endswith('.fq') or sys.argv[1].endswith('.fastq'):
        filetype = 'fq'
    if not filetype:
        filetype = sys.argv[2]
    continuously_write_output(seq_length(load_data(sys.argv[1]), sys.argv[2]), sys.stdout)
else:
    if os.path.exists(sys.argv[3]):
        os.remove(sys.argv[3])
    if sys.argv[1].endswith('.fa') or sys.argv[1].endswith('.fasta'):
        filteype = 'fa'
    elif sys.argv[1].endswith('.fq') or sys.argv[1].endswith('.fastq'):
        filetype = 'fq'
    if not filetype:
        filetype = sys.argv[2]
    continuously_write_output(seq_length(load_data(sys.argv[1]), sys.argv[2]), open(sys.argv[3], 'w'))
