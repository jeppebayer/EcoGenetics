#!/bin/bash

sacct -o AllocCPUs,ReqCPUS,MaxRSS,ReqMem,Elapsed,ElapsedRaw,Timelimit,JobID -n -P --delimiter=~ --units=M -j 11246544 | sort --field-separator=~ -n -k5 -r | less

# Remove large number of files
find people/Jeppe_Bayer/steps/03_sfs/out/ -maxdepth 1 -name "*.out" -print0 | xargs -0 rm

split()
{
    parallel \
    --pipepart \
    --block -1 \
    -k \
    -a /faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/Ocin_NYS-F/Ocin_NYS-F.pileup \
    grep \
    -F \
    "$region" \
    > /faststorage/project/EcoGenetics/people/Jeppe_Bayer/steps/temp/Ocin_NYS-F_"$num".mpileup
}