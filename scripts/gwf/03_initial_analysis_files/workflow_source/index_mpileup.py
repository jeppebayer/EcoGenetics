#!/bin/env python
import os, sys, re

def load_data(x: str):
    """Loads data either from standard input or from argument position 1.
    
    :param any x:
        Input data or '-'."""
    if x == '-':
        data = sys.stdin
    else:
        data = open(x, 'r')
    return data

# Unnecessary
# def reader(data):
#     for row in data:
#         yield row

def parse_mpileup(mpileup: str):
    """Parses :format:`mpileup` file yielding names of contigs and their respective start and end position
    
    ::

        yield 'contig\tstart\tend'

    :param str mpileup:
        Input :format:`mpileup` file."""
    yield 'contig\tstart\tend'
    contig = None
    line_num = 0
    rows = (row for row in mpileup) # Generator expression which yields each row in mpileup
    entries = (re.split(r'\s+', entry.rstrip(' '))[0] for entry in rows) # Generator expression which removes trailing spaces and splits rows on blank space characters and yields the first column entry of each row
    for entry in entries:
        line_num += 1
        if entry != contig:
            if contig:
                yield "{}\t{}\t{}".format(contig, str(start), str(line_num - 1))
            start = line_num
            contig = entry
    yield "{}\t{}\t{}".format(contig, str(start), str(line_num))
    
def continuously_write_output(input: str, output: str):
    """Continuously writes each iteration of input to output"""
    with output as outfile:
        for entry in input:
            outfile.write("{}\n".format(entry))
            outfile.flush()

if len(sys.argv) <= 1 or len(sys.argv) >= 4:
    sys.stdout.write('ERROR: {}: incorrect number of arguments\n'.format(os.path.basename(sys.argv[0])))
    exit(1)
elif len(sys.argv) == 2:
    continuously_write_output(parse_mpileup(load_data(sys.argv[1])), sys.stdout)
else:
    if os.path.exists(sys.argv[2]):
        os.remove(sys.argv[2])
    continuously_write_output(parse_mpileup(load_data(sys.argv[1])), open(sys.argv[2], 'w'))
