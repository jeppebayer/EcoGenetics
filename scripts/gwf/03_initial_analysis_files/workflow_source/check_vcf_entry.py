#!/bin/env python
from sys import stdin, stdout, argv 

if len(argv) > 1:
    in_ln = open(argv[1], 'r')
else:
    in_ln = stdin

count = 0
for ln in in_ln:
    ln = ln.strip()
    if not ln.startswith('#'):
        count += 1
    if count > 0 :
        break
stdout.write('{}\n'.format(count))

if len(argv) > 1:
    in_ln.close()