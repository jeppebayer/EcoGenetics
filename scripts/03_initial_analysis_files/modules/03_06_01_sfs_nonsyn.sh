#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal

length=${#$SLURM_ARRAY_TASK_ID}

if [ "$length" -lt 2 ]; then
    num="0$SLURM_ARRAY_TASK_ID"
else
    num="$SLURM_ARRAY_TASK_ID"
fi

exit 0