#!/bin/bash

srun --cpus-per-task=8 --mem=20g --time=06:00:00 --account=EcoGenetics --pty bash

java -jar -Xmx20g "$(readlink -f "$(dirname "$0")")"/Juicebox_*.jar
