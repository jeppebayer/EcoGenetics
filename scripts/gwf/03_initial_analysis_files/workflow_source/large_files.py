#!/bin/env python
import os, sys, re

def load_data(x):
    if x == '-':
        data = sys.stdin
    else:
        data = open(x, 'r')
    return data

def parse_mpileup(mpileup: str):
    contig = None
    line_num = 0
    for entry in mpileup:
        line_num += 1
        entry = re.split(r'\s+', entry.rstrip(' '))[0]
        if entry != contig:
            if contig:
                yield {'contig': contig, 'start': str(start), 'end': str(line_num - 1)}
            start = line_num
            contig = entry
    yield {'contig': contig, 'start': str(start), 'end': str(line_num)}

def write_dict_to_table(data, output: str):
    if os.path.exists(output):
        os.remove(output)
    header = None
    with open(output, 'w') as outfile:
        for entry in data:
            if not header:
                header = outfile.write('{}\n'.format('\t'.join(list(entry.keys()))))
            outfile.write('{}\n'.format('\t'.join(list(entry.values()))))

write_dict_to_table(parse_mpileup(load_data(sys.argv[1])), sys.argv[2])