#!/bin/bash
#SBATCH -c 4 -t 1-0

# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euo pipefail

printf "Running snakemake...\n"
snakemake \
    -j 4 \
    --software-deployment-method conda \
    --rerun-incomplete
printf "Run of snakemake complete.\n"
