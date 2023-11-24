#!/bin/env python
import os, sys, gzip

usage = "\nUsage: {} fa|fq [input file] [output file]\n".format(sys.argv[0])

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
        print('Unrecognized filetype: \'{}\''.format(filetype))
        exit(1)
    rows = (row for row in sequence_file)
    entries = (entry.rstrip() for entry in rows)
    if filetype == 'fa':
        seq_name = None
        length = 0
        for entry in entries:
            if entry.startswith('>'):
                if seq_name:
                    yield "{}\t{}".format(seq_name, length)
                    length = 0
                seq_name = entry.split(" ", 1)[0][1:]
            else:
                length += len(entry)
        yield "{}\t{}".format(seq_name, length)
    if filetype == 'fq':
        for entry in entries:
            if entry.startswith('@'):
                yield "{}\t{}".format(entry[1:], len(next(entries)))
                next(entries)
                next(entries)

def continuously_write_output(input: str, output: str):
    """Continuously writes each iteration of input to output"""
    with output as outfile:
        for entry in input:
            outfile.write("{}\n".format(entry))
            outfile.flush()

if len(sys.argv) <= 2 or len(sys.argv) >= 5:
    sys.stdout.write('{}\n'.format(usage))
    exit(1)
elif len(sys.argv) == 3:
    continuously_write_output(seq_length(load_data(sys.argv[1]), sys.argv[2]), sys.stdout)
else:
    if os.path.exists(sys.argv[3]):
        os.remove(sys.argv[3])
    continuously_write_output(seq_length(load_data(sys.argv[1]), sys.argv[2]), open(sys.argv[3], 'w'))
