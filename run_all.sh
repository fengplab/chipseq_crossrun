#!/usr/bin/env bash
# Run the cross-run summary end to end, or selected steps.
# Usage: ./run_all.sh config.yaml [first_step [last_step]]
#   steps: 0 setup | 1 inputs | 2 motif discovery (XSTREME with control-library sequences, or the pipeline's)
#          3 sequences | 4 motif families | 5 motif scans
#          6 repeat enrichment | 7 reads in repeats (optional) | 8 regions
#          9 R: motifs (Q1) | 10 R: repeats (Q2) | 11 R: known motifs (Q3) | 12 R: regions (Q4) | 13 R: metadata (Q5)
#          14 QC computations (mappability, insert sizes) | 15 R: QC | 16 report (all figures + data dictionary)
# Steps are idempotent: finished per-sample outputs are reused, delete them to force a rerun.
set -euo pipefail
CFG="$(cd "$(dirname "${1:?config.yaml}")" && pwd)/$(basename "$1")"
FIRST="${2:-0}"; LAST="${3:-16}"
ROOT="$(cd "$(dirname "$0")" && pwd)"
OUT=$(python3 -c "import yaml,os,sys; print(os.path.abspath(yaml.safe_load(open(sys.argv[1]))['outdir']))" "$CFG")
ENV="$OUT/work/config.env"

# --- R isolation ---------------------------------------------------------------------------------
# When R comes from the active conda env, stop it from also searching the user's personal library
# (~/R/<platform>-library/<version>). Packages there are built for a different R installation and can
# fail or segfault when conda's R loads them. Set CROSSRUN_ALLOW_USER_RLIB=1 to disable this.
RSCRIPT="$(command -v Rscript || true)"
if [[ -n "${CONDA_PREFIX:-}" && "$RSCRIPT" == "$CONDA_PREFIX"/* && "${CROSSRUN_ALLOW_USER_RLIB:-0}" != "1" ]]; then
    export R_LIBS_USER="$CONDA_PREFIX/.crossrun_no_user_rlib"   # non-existent path -> ignored by R
    export R_LIBS=""
fi

# --- R preflight: load every package before any slow step runs ----------------------------------
if (( LAST >= 4 )); then
    [[ -n "$RSCRIPT" ]] || { echo "[preflight] Rscript not found on PATH" >&2; exit 1; }
    if ! "$RSCRIPT" -e '
        pk <- c("data.table", "ggplot2", "yaml", "pheatmap", "UpSetR", "patchwork", "lme4", "lmerTest", "vegan")
        for (p in pk) suppressPackageStartupMessages(library(p, character.only = TRUE))
        cat("[preflight] R packages OK; libraries:", paste(.libPaths(), collapse = " "), "\n")' >&2; then
        echo "[preflight] R failed to load its packages (see above)." >&2
        echo "[preflight] Check which copy is being loaded:  Rscript -e 'print(.libPaths())'" >&2
        exit 1
    fi
fi

step() {
    local n="$1"; shift
    (( n < FIRST || n > LAST )) && return 0
    echo "================ step $n: $1 ================" >&2
    shift
    mkdir -p "$OUT/logs"
    "$@" 2>&1 | tee "$OUT/logs/step$(printf %02d "$n").log" >&2
    return "${PIPESTATUS[0]}"
}

step 0  "setup"                 python3 "$ROOT/bin/00_setup.py" "$CFG"
step 1  "prepare inputs"        bash "$ROOT/steps/01_prepare_inputs.sh" "$ENV"
step 2  "motif discovery"       bash "$ROOT/steps/02_motif_discovery.sh" "$ENV"
step 3  "scan sequences"        bash "$ROOT/steps/03_sequences.sh" "$ENV"
step 4  "motif families"        bash "$ROOT/steps/04_motif_families.sh" "$ENV"
step 5  "motif scans"           bash "$ROOT/steps/05_motif_scan.sh" "$ENV"
step 6  "repeat enrichment"     bash "$ROOT/steps/06_repeat_enrichment.sh" "$ENV"
step 7  "reads in repeats"      bash "$ROOT/steps/07_reads_in_repeats.sh" "$ENV"
step 8  "regions"               bash "$ROOT/steps/08_regions.sh" "$ENV"
step 9  "Q1 motifs"             Rscript "$ROOT/R/10_motifs.R" "$CFG"
step 10 "Q2 repeats"            Rscript "$ROOT/R/20_repeats.R" "$CFG"
step 11 "Q3 known motifs"       Rscript "$ROOT/R/30_known_motifs.R" "$CFG"
step 12 "Q4 regions"            Rscript "$ROOT/R/40_regions.R" "$CFG"
step 13 "Q5 metadata"           Rscript "$ROOT/R/50_metadata.R" "$CFG"
step 14 "QC computations"       bash "$ROOT/steps/14_qc.sh" "$ENV"
step 15 "QC analysis"           Rscript "$ROOT/R/60_qc.R" "$CFG"
step 16 "report"                python3 "$ROOT/bin/build_report.py" "$CFG"
echo "done: tables in $OUT/tables, figures in $OUT/figures, report in $OUT/report/crossrun_report.html" >&2
