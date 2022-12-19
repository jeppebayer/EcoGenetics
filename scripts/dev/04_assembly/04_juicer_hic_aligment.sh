#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

bash /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/juicer/scripts/juicer.sh \
-d \
-D \
-p chrom.sizes \
-s none \
-z [fasta] \
-q short \
-Q 12:00:00 \
-l normal \
-L 24:00:00 \
-t 36 \
> [species]_juicer.log

exit 0