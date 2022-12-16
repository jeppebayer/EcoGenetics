#!/bin/bash

sacct -o AllocCPUs,ReqCPUS,MaxRSS,ReqMem,Elapsed,ElapsedRaw,Timelimit,JobID -n -P --delimiter=~ --units=M -j 11246544 | sort --field-separator=~ -n -k5 -r | less

# Remove large number of files
find people/Jeppe_Bayer/steps/03_sfs/out/ -maxdepth 1 -name "*.out" -print0 | xargs -0 rm