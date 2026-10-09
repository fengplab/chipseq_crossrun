#!/usr/bin/env bash
# Smoke test on synthetic data with stub fimo/tomtom/sea (exact-consensus matching; NOT the real tools).
# Needs bedtools, python3 (pyyaml), R with the packages listed in environment.yml.
set -euo pipefail
T="$(cd "$(dirname "$0")" && pwd)"; ROOT="$T/.."
D="${1:-$T/test_data}"
python3 "$T/make_test_data.py" "$D"
sed -e "s#manifest: manifest.tsv#manifest: $D/manifest.tsv#" -e "s#outdir: crossrun_out#outdir: $D/out#" \
    -e "s#genome_build: T2T-CHM13v2.0#genome_build: test#" -e "s#/path/to/genome.fa.sizes#$D/genome.sizes#" \
    -e "s#/path/to/genome.fa#$D/genome.fa#" -e "s#include_bed: null#include_bed: $D/include.bed#" \
    -e "s#/path/to/repeatmasker.bb#$D/rmsk.bed#" -e "s#/path/to/censat.bb#$D/censat.bed#" \
    -e "s#genes_bed: null #genes_bed: $D/genes.bed#" -e "s#  - /path/to/JASPAR2024_CORE_vertebrates_non-redundant.meme#  - $D/known.meme#" \
    -e "/H12CORE_meme_format/d" -e "s#n_shuffles: 200#n_shuffles: 20#" "$ROOT/config.yaml" > "$D/config.yaml"
PATH="$T/fakebin:$PATH" "$ROOT/run_all.sh" "$D/config.yaml"

# Checks on the end-to-end run. The stub FIMO renames sequences like MEME Suite >= 5.5.4 does unless --no-pgc is used
# or headers carry no coordinates, so these fail if motif sites stop reaching the peaks.
python3 - "$D/out" <<'PY'
import csv, gzip, glob, sys
o = sys.argv[1]
rd = lambda f: list(csv.DictReader(open(f"{o}/tables/{f}.tsv"), delimiter="\t"))
sites = sum(1 for f in glob.glob(f"{o}/work/scan/*.peaks.sites.bed.gz") for _ in gzip.open(f, "rt"))
assert sites > 0, "no motif sites recorded in peaks (step 05)"
t0 = rd("table00_peak_summary")
assert t0 and all(int(r["n_peaks"]) > 0 for r in t0), "peak overview: a sample has no peaks in table00"
assert any(r["fold_enrichment_median"] not in ("NA", "") for r in t0), "peak overview: no fold enrichment parsed"
t0b = rd("table00b_peaks_by_chromosome")
for r in t0:
    n = sum(int(x["n_peaks"]) for x in t0b if x["sample_id"] == r["sample_id"])
    assert n == int(r["n_peaks"]), f"peak overview: per-chromosome counts do not add up for {r['sample_id']}"
for f in ("pk01_peaks_per_sample", "pk02_peaks_by_chromosome", "pk03_peak_width", "pk04_fold_enrichment"):
    assert glob.glob(f"{o}/figures/{f}.png"), f"missing figure {f}"
t9 = rd("table09_motif_classification")
assert any(int(r["n_motif_peaks"] or 0) > 0 for r in t9), "no motif family has a site in any modelled peak (step 16/17)"
t12 = rd("table12_timepoint_contrasts")
bad = [r for r in t12 if r["contrast_se"] not in ("NA", "") and float(r["contrast_se"]) == 0 and r["q_value"] not in ("NA", "")]
assert not bad, f"{len(bad)} time-point contrasts with zero standard error were tested"
for t in ("table06_metadata_effects", "table07_motif_detection_vs_metadata", "table07b_permanova", "table06c_timecourse_trends"):
    assert not glob.glob(f"{o}/tables/{t}.tsv"), f"{t} should no longer be written"
print(f"[test] end-to-end: {sites} motif sites in peaks; {sum(1 for r in t9 if int(r['n_motif_peaks'] or 0) > 0)} families with motif peaks")
PY

# Regression models (R/70_models.R) on simulated windows with planted effects: checks the classifications and the
# planted covariate effects are recovered. Needs only R.
M="${D}_models"
python3 "$T/make_model_test_data.py" "$M" 600
Rscript "$ROOT/R/70_models.R" "$M/config.yaml"
python3 - "$M/out/tables" <<'PY'
import csv, sys
t = sys.argv[1]
rd = lambda f: list(csv.DictReader(open(f"{t}/{f}.tsv"), delimiter="\t"))
cls = {r["family_id"]: r["classification"] for r in rd("table09_motif_classification")}
want = {"F001": "universal", "F002": "repeat-driven", "F003": "repeat-specific binding", "F004": "confined to repeats",
        "F006": "universal", "F007": "universal"}
bad = {f: (cls.get(f), w) for f, w in want.items() if cls.get(f) != w}
assert not bad, f"classification mismatch (got, expected): {bad}"
assert "F005" not in cls, "F005 should fail the family cutoff"
sig = {(r["feature_type"], r["response_id"], r["covariate"]) for r in rd("table11_covariate_tests") if r["q_value"] != "NA" and float(r["q_value"]) < 0.05}
for k in [("density", "F006", "virus"), ("density", "F007", "time_bin"), ("peak_share", "censat:hsat3", "time_bin"), ("reads_class", "censat:hsat3", "time_bin")]:
    assert k in sig, f"planted effect not detected: {k}"
cv = {(r["family_id"], r["covariate"]) for r in rd("table09c_peak_model_covariates") if float(r["q_value"]) < 0.05}
assert cv == {("F006", "virus"), ("F007", "time_bin")}, f"covariate x peak effects: {cv}"
print("[test] regression models: planted classifications and covariate effects recovered;", len(sig), "sample-level tests with q < 0.05")
PY
