#!/usr/bin/env python3
"""
Simulated inputs for R/70_models.R with known answers (no genome, no external tools).

Writes <dir>/config.yaml and <dir>/out/work/... exactly as steps 0-16 would, for 24 samples:
  virus A, B, C x cell type X, Y (virus C only in cell type Y), three time points (2, 8, 24) per run, two runs per
  virus x cell type; each run is one time-course series.

Planted motif families (expected classification by the peak-level model):
  F001 universal                 5% of background windows, 25% of peaks, in every repeat class
  F002 repeat-driven             2% outside repeats, 50% inside active HOR, peaks and background alike
  F003 repeat-specific binding   5% everywhere, except 40% of peaks inside SINE/Alu
  F004 confined to repeats       98% of HSat3 windows, nowhere else; long tandem arrays
  F005 no association            8% everywhere
  F006 universal, virus effect   5% background; peaks 10% (virus A, C) or 35% (virus B)
  F007 universal, time effect    5% background; peaks 10% / 20% / 35% at 2 / 8 / 24 h
Planted sample-level effects: share of peaks in HSat3 and IP reads in HSat3 rise with time; virus B has more
F006 per kb; F007 density rises with time. Everything else is noise.

Usage: make_model_test_data.py <dir> [peaks_per_sample]
"""
import gzip
import os
import sys

import numpy as np

D = os.path.abspath(sys.argv[1]); NPK = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
rng = np.random.default_rng(7)
W = f"{D}/out/work"
for d in ("models/win", "motifs", "scan", "reads", "peaks"):
    os.makedirs(f"{W}/{d}", exist_ok=True)
os.makedirs(f"{D}/out/tables", exist_ok=True)

samples = []
run = 0
for virus, cell in [("A", "X"), ("A", "Y"), ("B", "X"), ("B", "Y"), ("C", "Y")]:
    for rep in range(2):
        run += 1
        for t in (2, 8, 24):
            samples.append(dict(sample_id=f"{virus}_{cell}_r{run}_{t}h_REP1", run_id=f"run{run:02d}", series_id=f"{virus}_{cell}_r{run}",
                                virus=virus, cell_type=cell, hpi=t, cgas_construct="WT"))
samples = samples[:24] + samples[24:]
with open(f"{W}/samples.tsv", "w") as fh:
    cols = ["sample_id", "run_id", "series_id", "virus", "cell_type", "hpi", "cgas_construct", "peaks_file", "bam_file", "control_bam"]
    fh.write("\t".join(cols) + "\n")
    for s in samples:
        s.update(peaks_file="x", bam_file=f"/bam/{s['sample_id']}.bam", control_bam=f"/bam/INPUT_{s['run_id']}.bam")
        fh.write("\t".join(str(s[c]) for c in cols) + "\n")

CLASSES = ["none", "censat:hsat3", "censat:active_hor", "repeatmasker:SINE/Alu"]
FAMS = {"F001": "ACGTGACT", "F002": "TTCGTTGGAAACGGGA", "F003": "GGCCGGGCG", "F004": "TTCCATTCCA", "F005": "ATGCATAA",
        "F006": "CCGGAAGT", "F007": "GATTACAG"}


def p_present(fam, peak, cls, s):
    tb = {2: 0, 8: 1, 24: 2}[s["hpi"]]
    if fam == "F001": return 0.25 if peak else 0.05
    if fam == "F002": return 0.50 if cls == "censat:active_hor" else 0.02
    if fam == "F003": return 0.40 if (peak and cls == "repeatmasker:SINE/Alu") else 0.05
    if fam == "F004": return 0.98 if cls == "censat:hsat3" else 0.0
    if fam == "F005": return 0.08
    if fam == "F006": return (0.35 if s["virus"] == "B" else 0.10) if peak else 0.05
    if fam == "F007": return (0.10, 0.20, 0.35)[tb] if peak else 0.05


lib = open(f"{W}/models/libsize.tsv", "w"); lib.write("sample_id\tip_total_reads\tctl_total_reads\n")
rl = open(f"{W}/reads/reads_in_labels.tsv", "w"); rl.write("bam\tfeature_set\tlabel\treads\ttotal_reads\n")
done_ctl = set()
for s in samples:
    sid = s["sample_id"]; tb = {2: 0, 8: 1, 24: 2}[s["hpi"]]
    hs = (0.06, 0.12, 0.22)[tb]                      # share of peaks in HSat3 rises with time
    pk_p = [1 - hs - 0.18, hs, 0.08, 0.10]; bg_p = [0.80, 0.05, 0.05, 0.10]
    rows, sites = [], []
    slope = rng.normal(0, 0.25)                      # sample-specific shift of the peak effect (logit scale)
    for i in range(NPK):
        for draw in range(4):                        # 0 = the peak, 1-3 = its background windows
            peak = draw == 0
            wid = f"{sid}|peak_{i}" + ("" if peak else f"#bg{draw}")
            cls = CLASSES[rng.choice(4, p=pk_p if peak else bg_p)]
            width = int(rng.integers(200, 600))
            ip = rng.poisson(30 * np.exp(rng.normal(0, 0.5))) if peak else rng.poisson(5)
            ctl = rng.poisson(8)
            cs, rm = ("none", "none")
            if cls.startswith("censat:"): cs = cls.split(":", 1)[1]
            if cls.startswith("repeatmasker:"): rm = cls.split(":", 1)[1]
            rows.append(f"{sid}\t{wid}\t{'peak' if peak else 'background'}\t{draw}\t{sid}|peak_{i}\t1\tchr1\t{i*1000}\t{i*1000+width}\t{width}\t"
                        f"{rng.normal(0.42, 0.05):.4f}\t{'NA' if not peak else '5'}\t{'R%d' % (i % (NPK // 2)) if peak else 'NA'}\t"
                        f"{cs}\t{1.0 if cs != 'none' else 0.0}\t{rm}\t{1.0 if rm != 'none' else 0.0}\t{ip}\t{ctl}\n")
            for fam in FAMS:
                p = p_present(fam, peak, cls, s)
                if 0 < p < 1 and peak and fam in ("F001", "F006", "F007"):
                    p = 1 / (1 + np.exp(-(np.log(p / (1 - p)) + slope)))
                if rng.random() < p:
                    n = 1 + rng.poisson(0.2) if fam != "F004" else int(rng.integers(15, 60))
                    tand = n if fam == "F004" else (n if n >= 3 and rng.random() < 0.3 else 0)
                    sites.append(f"{sid}\t{wid}\t{fam}\t{n}\t{tand}\n")
    with gzip.open(f"{W}/models/win/{sid}.windows.tsv.gz", "wt") as fh:
        fh.write("sample_id\twindow_id\tset\tbg_draw\tpeak_id\tin_model\tchrom\tstart\tend\twidth\tgc\tneglog10_q\tregion_id\t"
                 "censat_class\tcensat_frac\trm_class\trm_frac\tip_reads\tctl_reads\n")
        fh.writelines(rows)
    with gzip.open(f"{W}/models/win/{sid}.sites.tsv.gz", "wt") as fh:
        fh.write("sample_id\twindow_id\tfamily_id\tn_sites\tn_tandem_sites\n"); fh.writelines(sites)
    lib.write(f"{sid}\t1000000\t1000000\n")
    for fs, lab, frac in [("censat", "hsat3", (0.02, 0.04, 0.08)[tb]), ("censat", "active_hor", 0.03), ("repeatmasker", "SINE/Alu", 0.10)]:
        rl.write(f"{s['bam_file']}\t{fs}\t{lab}\t{rng.binomial(1000000, frac * np.exp(rng.normal(0, 0.15)))}\t1000000\n")
    if s["control_bam"] not in done_ctl:
        done_ctl.add(s["control_bam"])
        for fs, lab, frac in [("censat", "hsat3", 0.02), ("censat", "active_hor", 0.03), ("repeatmasker", "SINE/Alu", 0.10)]:
            rl.write(f"{s['control_bam']}\t{fs}\t{lab}\t{rng.binomial(1000000, frac * np.exp(rng.normal(0, 0.05)))}\t1000000\n")
    with open(f"{W}/scan/{sid}.family_sea.tsv", "w") as fh:
        fh.write("RANK\tDB\tID\tALT_ID\tCONSENSUS\tTP\tTP%\tFP\tFP%\tENR_RATIO\tSCORE_THR\tPVALUE\tLOG_PVALUE\tEVALUE\tLOG_EVALUE\tQVALUE\tLOG_QVALUE\n")
        for k, fam in enumerate(FAMS):
            p = 0.4 if fam == "F005" else 10.0 ** -rng.uniform(3, 12)
            fh.write(f"{k+1}\tdb\t{fam}\tx\t{FAMS[fam]}\t1\t1\t1\t1\t1\t1\t{p:.3g}\t{np.log(p):.2f}\t{p:.3g}\t0\t{p:.3g}\t0\n")
        fh.write("\n# SEA (Simple Enrichment Analysis): stub\n")
lib.close(); rl.close()
with open(f"{W}/motifs/families.tsv", "w") as fh:
    fh.write("motif_uid\tsample_id\tmotif_id\tconsensus\twidth\tsource\tsea_pvalue\tevalue\tis_representative\tfamily_id\tfamily_n_samples\tfamily_size\n")
    for fam, c in FAMS.items():
        for s in samples[:3]:
            fh.write(f"{s['sample_id']}__{c}\t{s['sample_id']}\t{c}\t{c}\t{len(c)}\tSTREME\t1e-8\t1e-6\tFALSE\t{fam}\t3\t3\n")
with open(f"{W}/motifs/family_rep_features.tsv", "w") as fh:
    fh.write("family_id\trep_motif_uid\tpwm_consensus\twidth\tinfo_content\tself_period_shift\tself_period_corr\ttrimer_diversity\n")
    for fam, c in FAMS.items():
        fh.write(f"{fam}\tx__{c}\t{c}\t{len(c)}\t12\t1\t0.1\t1\n")
with open(f"{D}/config.yaml", "w") as fh:
    fh.write(f"""outdir: {D}/out
threads: 8
batch_var: run_id
series_var: series_id
time_var: hpi
categorical_covariates: [virus, cell_type, cgas_construct]
numeric_covariates: [hpi]
covariate_labels: {{hpi: "Hours post-infection", cell_type: "Cell type", cgas_construct: "cGAS construct", virus: "Virus"}}
models:
  family_stat: sea_pvalue
  family_cutoff: 0.05
  time_bins: [4, 12]
  n_permutations: 199
  n_bootstrap: 100
""")
print(f"[make_model_test_data] {len(samples)} samples x {NPK} peaks -> {D}")
