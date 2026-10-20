#!/bin/bash
#SBATCH --job-name=compare_pies
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=compare_pies_%j.out
# One-off pie comparisons between groups of samples (extra/compare_pies.R). Runs setup and input preparation
# (run_all.sh steps 0-1), then draws every comparison in the TSV. IP-coverage pies are added when the main run's
# step 9 output (work/contigs/contig_reads.tsv) exists. Submit from the bundle directory with the crossrun env active:
#   sbatch extra/compare_pies.sbatch /path/to/config.yaml extra/comparisons_cgas.tsv [--equal-weight] [--peaks-only]
# Output: <outdir>/extra/pies/compare/
set -euo pipefail
CFG="${1:?config.yaml}"; CMP="${2:?comparisons.tsv}"; shift 2
bash run_all.sh "$CFG" 0 1
Rscript extra/compare_pies.R "$CFG" "$CMP" "$@"
