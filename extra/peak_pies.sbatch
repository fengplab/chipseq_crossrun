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
#   sbatch extra/peak_pies.sbatch /path/to/config.yaml [virus,cell_type]
# The optional second argument pools samples by those manifest columns (one figure per column, one pie per level),
# in addition to extra_plots.pie_group_by in the config: extra/pies/peaks_by_<column>_{ring,noring}.{png,tsv}
# Output: <outdir>/extra/pies/peaks_<sample>_{ring,noring}.{png,tsv}, pooled_peaks_{ring,noring}.{png,tsv},
#         peaks_per_chromosome.tsv (chromosome x sample counts)
set -euo pipefail
CFG="${1:?config.yaml}"
bash run_all.sh "$CFG" 0 1
GROUP="${2:-}"
Rscript R/plot_pies.R "$CFG" --peaks-only ${GROUP:+--group-by "$GROUP"}
