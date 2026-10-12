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
# viral genomes and chrM: kept out of the peak analyses, counted in the read breakdowns
gs = [l.split("\t")[0] for l in open(f"{o}/work/annot/genome.sizes")]
assert "chrM" not in gs and not any(x in gs for x in ("HSV-1", "KSHV", "ADV")), f"genome.sizes should hold host nuclear chromosomes only: {gs}"
ga = [l.split("\t")[0] for l in open(f"{o}/work/annot/genome_all.sizes")]
assert all(x in ga for x in ("chrM", "HSV-1", "KSHV", "ADV")), f"genome_all.sizes misses contigs: {ga}"
c0 = rd("table00c_contig_reads")
def reads(sample_prefix, lib, contig, col="reads"):
    return sum(int(r[col]) for r in c0 if r["sample_id"].startswith(sample_prefix) and r["library"] == lib and r["contig"] == contig)
assert reads("HSV_", "IP", "HSV-1") > 0, "no HSV-1 reads counted in HSV-1 IP libraries"
assert reads("HSV_", "IP", "KSHV", "reads_all") > 0 and reads("HSV_", "IP", "KSHV") == 0, "low-MAPQ KSHV reads should be counted before but not after the MAPQ filter"
assert reads("HSV_", "IP", "HSV-1") < reads("HSV_", "IP", "HSV-1", "reads_all"), "the viral MAPQ filter removed nothing from HSV-1"
assert reads("VACV_", "IP", "ADV") > 0 and reads("HSV_", "IP", "ADV") == 0, "adenovirus reads expected only in HEK293T libraries"
assert reads("HSV_", "IP", "chrM") > 0, "no chrM reads counted"
q7 = rd("qc07_viral_mapq_summary")
hq = [r for r in q7 if r["contig"] == "HSV-1" and r["library"] == "IP" and int(r["reads_before_filter"]) > 0]
assert hq and all(float(r["pct_passing"]) < 100 for r in hq), "qc07: HSV-1 rows should show reads below the MAPQ cutoff"
assert glob.glob(f"{o}/figures/pk05_viral_reads.png"), "missing Figure P5"
assert not glob.glob(f"{o}/figures/qc03_insert_*"), "QC3 figures should no longer be in figures/"
assert glob.glob(f"{o}/extra/insert_sizes/*.png"), "insert-size figures missing from extra/insert_sizes"
s1 = "HSV_THP1_WT_24h_REP1"
for f in ("chr1_raw", "chr1_cpm", "chr1_genome", "chr1_input", "chrM_raw", "HSV-1_raw", "HSV-1_input", "ADV_genome"):
    assert glob.glob(f"{o}/extra/coverage/{s1}/{f}.png"), f"missing coverage plot {s1}/{f}"
n_cov = len(glob.glob(f"{o}/extra/coverage/*/*.png"))
assert n_cov == 10 * 8 * 4, f"expected 320 coverage plots (10 samples x 8 contigs x 4 normalisations), found {n_cov}"
for f in (f"coverage_{s1}_with_viruses_ring", f"coverage_{s1}_with_viruses_noring", f"coverage_{s1}_host_ring", f"coverage_{s1}_host_noring",
          f"peaks_{s1}_ring", f"peaks_{s1}_noring", "pooled_coverage_rpkm_with_viruses", "pooled_coverage_rpkm_host", "pooled_peaks_ring", "pooled_peaks_noring"):
    assert glob.glob(f"{o}/extra/pies/{f}.png") and glob.glob(f"{o}/extra/pies/{f}.tsv"), f"missing pie {f}"
assert not glob.glob(f"{o}/extra/pies/pooled_coverage_rpkm_*ring*"), "the pooled RPKM pie has no ring version"
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
