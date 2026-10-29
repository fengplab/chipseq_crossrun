#!/bin/bash
#SBATCH --job-name=level1
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=level1_%A_%a.out
# Binding models, Level 1: one array task per line of the task file (sample_id <TAB> arm), written by step 18.
# Submit from the bundle directory with the crossrun environment active (step 18 prints the exact command):
#   sbatch --array=1-<N> models/level1.sbatch /path/to/config.yaml <outdir>/work/models/level1_tasks_todo.tsv
# Size --cpus-per-task / --mem from tests/bench_level1.sbatch at your bin size (memory grows with bins x families);
# for the big-memory nodes use --partition=bigmem on the sbatch command line. Set binding_models.level1_threads to
# the CPUs given here.
set -euo pipefail
CFG="${1:?config.yaml}"; TASKS="${2:?task file}"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID:?run as an array job}p" "$TASKS")
[[ -n "$LINE" ]] || { echo "no task $SLURM_ARRAY_TASK_ID in $TASKS" >&2; exit 1; }
IFS=$'\t' read -r SAMPLE ARM <<< "$LINE"
echo "task $SLURM_ARRAY_TASK_ID: $SAMPLE $ARM" >&2
Rscript R/70_level1.R "$CFG" "$SAMPLE" "$ARM"
