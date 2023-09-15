#!/bin/env python
from sys import stdin, stdout, argv 

def check_vcf(in_ln: str):
    """Check if :format:`VCF` file contains more than header.
    
    This funtions simply check if there are any non-whitespace lines after the file header.
    As soon as one if found it returns 1 (int), if non is found returns 0 (int).
    
    :param str in_ln:
        :format:`VCF` content."""
    count = 0
    for ln in in_ln:
        ln = ln.strip()
        if not ln.startswith('#'):
            count += 1
        if count > 0 :
            break
    return count

# Run check_vcf on either stdin or file argument
if len(argv) > 1:
    in_ln = open(argv[1], 'r')
else:
    in_ln = stdin

stdout.write('{}\n'.format(check_vcf(in_ln)))

# Remember to close file
if len(argv) > 1:
    in_ln.close()