#!/bin/env python
import os, sys, re

def load_data(x):
    """Loads data either from standard input or from argument position 1.
    
    :param any x:
        Input data or '-'."""
    if x == '-':
        data = sys.stdin
    else:
        data = open(x, 'r')
    return data

def parse_mpileup(mpileup: str):
    """Parses :format:`mpileup` yielding dictionaries containing contig name and corresponding start and end line.
    
    ::
    
        yield {'contig': str, 'start': int, 'end': int}
    
    :param str mpileup:
        Input :format:`mpileup` file."""
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

def listdict_to_table(data: list, output: str):
    """Writes list of dictionaries to :format:`text` file as table.
    
    :param list data:
        List of dictionaries.
    :param str output:
        Output file path/name."""
    if os.path.exists(output):
        os.remove(output)
    header = None
    with open(output, 'w') as outfile:
        for entry in data:
            if not header:
                header = outfile.write('{}\n'.format('\t'.join(list(entry.keys()))))
            outfile.write('{}\n'.format('\t'.join(list(entry.values()))))

listdict_to_table(parse_mpileup(load_data(sys.argv[1])), sys.argv[2])