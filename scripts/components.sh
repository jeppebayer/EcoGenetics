#!/bin/bash

# ----------------- Functions --------------------------------------------

# Function to set dependency on job according macthing the number of 
# relevant executed jobs. The function take a list of job numbers
# separated by spaces as input in a single argument.
dependency()
{
    depend_on=
    for id in $1; do
        if [ -n "${joblist[$id]}" ]; then
            if [ -z "$depend_on" ]; then
                depend_on="${joblist[$id]}"
            else
                depend_on="$depend_on:${joblist[$id]}"
            fi
        fi
    done
    [ -z "$depend_on" ] && depend_on="-1"
    echo "$depend_on"
}

# Function used to adjust time for jobs. Takes three argument as input.
# First; an adjustment value, which sets the time in relative terms
# to some standard, second; a base time in minutes, third; a static time
# in minutes, to serve as a lower bound on time.
timer()
{
    adjbase=$(awk -v adjustment="$1" -v base="$2" 'BEGIN { print int( ( base * adjustment ) * 1.10 ) }')
    if [ "$adjbase" -lt "$3" ]; then 
        awk -v adjbase="$adjbase" -v static="$3" 'BEGIN { print int( adjbase + static ) }'
    else
        echo "$adjbase"
    fi
}

# ----------------- Jobs -------------------------------------------------

# Holds IDs of executed jobs in array indexed acording to job ID
joblist=()

job1()
{
    joblist[1]=$(sbatch \
                --parsable \
                --time="$(timer 1 2 3)" \
                --mem=2G \
                --cpus-per-task=1 \
                --output=%j.out \
                script1.sh)
}

job2()
{
    joblist[2]=$(sbatch \
                --parsable \
                --time="$(timer 1 2 3)" \
                --mem=2G \
                --cpus-per-task=1 \
                --dependency=afterany:"$(dependency "1")" \
                --output=%j.out \
                script2.sh)
}

# ----------------- Script Queue -----------------------------------------

# Job switch
if [  ]; then
    job1
fi

if [  ]; then
    job2
fi

exit 0