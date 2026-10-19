#!/bin/bash
#SBATCH --job-name=peak_pies
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=peak_pies_%j.out
# Peaks per host chromosome as pie charts, and nothing else: setup and input preparation (run_all.sh steps 0-1),
# then R/plot_pies.R --peaks-only. chrM and the viral genomes are not peak-calling contigs: step 1 leaves their
# peaks out (and logs a warning if the pipeline called any), so the pies cover host nuclear chromosomes only.
# Submit from the bundle directory with the crossrun conda env active:
#   sbatch extra/peak_pies.sbatch /path/to/config.yaml
# Output: <outdir>/extra/pies/peaks_<sample>_{ring,noring}.{png,tsv}, pooled_peaks_{ring,noring}.{png,tsv},
#         peaks_per_chromosome.tsv (chromosome x sample counts)
set -euo pipefail
CFG="${1:?config.yaml}"
bash run_all.sh "$CFG" 0 1
Rscript R/plot_pies.R "$CFG" --peaks-only
