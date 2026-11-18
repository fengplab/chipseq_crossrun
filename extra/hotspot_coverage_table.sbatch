#!/bin/bash
#SBATCH --job-name=hotspot_cov
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=hotspot_cov_%j.out
# Supplemental table: hotspot bins with read coverage, after the whole crossrun has finished.
# Submit from the chipseq_crossrun bundle directory on BeeGFS, with the crossrun conda environment active
# (python3 and samtools come from it):
#   sbatch extra/hotspot_coverage_table.sbatch --crossrun <outdir> [options of extra/hotspot_coverage_table.py]
# Reads the bins and samples from <outdir> and the BAMs listed in <outdir>/work/samples.tsv; writes
# <outdir>/tables/supplemental_hotspot_bins.{tsv,per_library.tsv,columns.tsv,run_info.tsv}.
# Memory is small (one bin list per library); the time goes on one counting pass over each BAM.
set -euo pipefail
B="${CROSSRUN_BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ -r "$B/extra/hotspot_coverage_table.py" ] || {
    echo "cannot find extra/hotspot_coverage_table.py under $B: submit from the bundle directory or set CROSSRUN_BUNDLE" >&2
    exit 1; }
command -v samtools > /dev/null || { echo "samtools not on PATH: activate the crossrun environment before sbatch" >&2; exit 1; }
python3 "$B/extra/hotspot_coverage_table.py" --threads "${SLURM_CPUS_PER_TASK:-4}" "$@"
