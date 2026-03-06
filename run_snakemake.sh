#!/bin/bash
# run_snakemake.sh
# -----------------
# Launches the Snakemake pipeline with 4 cores by default
snakemake --cores 4 "$@"
