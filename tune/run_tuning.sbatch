#!/bin/bash
#SBATCH --job-name=bin_tuning
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=600G
#SBATCH --output=bin_tuning_%j.out
# Bin-size tuning on bigmem (753 GB, 16 cores per node). Submit from the bundle directory with the crossrun
# conda env active:   sbatch tune/run_tuning.sbatch /path/to/config.yaml
# Set threads: 16 in the config so the per-library signal runs 16 at a time.
set -euo pipefail
bash tune/run_tuning.sh "${1:?config.yaml}"
