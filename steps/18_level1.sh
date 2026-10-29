#!/usr/bin/env bash
# Step 18: binding models, Level 1 (one fit per modelled sample and read version; R/70_level1.R).
#   work/models/level1_tasks.tsv            one line per task: sample_id <TAB> arm
#   work/models/level1/<sample>__<arm>.*    coef.tsv, valid.tsv, info.tsv (written by each task)
# binding_models.level1_mode
#   inline : the tasks run here, one after another (small data, tests)
#   slurm  : when tasks are missing or out of date, this step writes the task list, prints the sbatch command for
#            models/level1.sbatch and stops; rerun run_all.sh from step 18 once the array job has finished.
# A task is up to date when its info.tsv says done/skipped and is newer than the bin table, the design matrix and the config.
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[[ "${MODELS_ENABLED:-1}" == "1" ]] || { log "binding models disabled (binding_models.enabled: false); skipping"; exit 0; }
M="$WORK/models"; L="$M/level1"; mkdir -p "$L"
B="$WORK/bins/w${BINS_SIZE:?no bin size; set bins.size or bins.tuning_record and rerun step 0}"
[[ -s "$B/libraries.tsv" ]] || { log "ERROR: $B/libraries.tsv missing; run step 16 first"; exit 1; }
[[ -s "$M/design_matrix.tsv" ]] || { log "ERROR: $M/design_matrix.tsv missing; run step 17 first"; exit 1; }
read -r MODE ARMS < <(python3 - "$CONFIG_YAML" <<'PY'
import sys, yaml
bm = (yaml.safe_load(open(sys.argv[1])) or {}).get("binding_models") or {}
arms = bm.get("arms") or ["all", "uniq"]
print(bm.get("level1_mode", "slurm"), ",".join(arms))
PY
)
python3 - "$M/design_matrix.tsv" "$B/libraries.tsv" "$ARMS" > "$M/level1_tasks.tsv" <<'PY'
import csv, sys
dm = [r["sample_id"] for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")]
libs = list(csv.DictReader(open(sys.argv[2]), delimiter="\t"))
have = {(r["sample_id"], r["library"]) for r in libs}
miss = [s for s in dm if (s, "IP") not in have or (s, "control") not in have]
if miss: sys.exit("samples of the design matrix without IP or control libraries in the bin table: " + ", ".join(miss))
for s in dm:
    for a in sys.argv[3].split(","): print(f"{s}\t{a}")
PY
N=$(wc -l < "$M/level1_tasks.tsv")
stale=()
while IFS=$'\t' read -r s a; do
    f="$L/${s}__${a}.info.tsv"
    if [[ ! -s "$f" ]] || ! grep -qE $'(^|\t)(done|skipped[^\t]*)(\t|$)' "$f" || [[ "$B/manifest.tsv" -nt "$f" || "$M/design_matrix.tsv" -nt "$f" || "$CONFIG_YAML" -nt "$f" ]]; then
        stale+=("$s"$'\t'"$a")
    fi
done < "$M/level1_tasks.tsv"
log "Level 1: $N tasks (samples x read versions $ARMS); ${#stale[@]} to run"
(( ${#stale[@]} == 0 )) && exit 0
if [[ "$MODE" == "inline" ]]; then
    for t in "${stale[@]}"; do IFS=$'\t' read -r s a <<< "$t"; Rscript "$ROOT/R/70_level1.R" "$CONFIG_YAML" "$s" "$a"; done
    exit 0
fi
printf "%s\n" "${stale[@]}" > "$M/level1_tasks_todo.tsv"
log "Level 1 runs as a Slurm array job. From the bundle directory, with the crossrun environment active:"
log "    sbatch --array=1-${#stale[@]} models/level1.sbatch $CONFIG_YAML $M/level1_tasks_todo.tsv"
log "then rerun:  ./run_all.sh $CONFIG_YAML 18"
exit 2
