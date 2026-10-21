#!/usr/bin/env python3
"""Sequence logos of the top motif families of a crossrun output, styled like the logo panel of the cGAS figure
(information content in bits, 0-2; A green, C blue, G yellow, T red; a numbered pink circle beside each logo).

Families are ranked from tables/table02_motif_prevalence.tsv (one row per family x sample); each logo is the family's
representative motif from work/motifs/family_reps.meme.

  --rank-by detected    number of samples in which the family is detected in peaks (column 'detected'),
                        ties broken by the median % of peaks with a site (default)
  --rank-by prevalence  median % of peaks with a site ('pct_peaks')
  --rank-by denovo      number of samples in which a member was discovered de novo ('discovered_denovo')
  --rank-by sea         smallest SEA q-value across samples ('sea_qvalue')
  --samples a,b         rank on these samples only (default: all samples in the table)
  --families F001,F007  draw exactly these families, in this order (no ranking)

Usage:
  python3 bin/plot_family_logos.py --crossrun <outdir> [--top 5] [--rank-by detected] [--samples ...]
                                   [--ncol 1] [--positions 30] [--width 7.5] [--row-height 0.9] [--dpi 300] [--out family_logos.png]
Writes the figure (.png, .pdf or .svg) and <out>.tsv with the families drawn, their rank and the ranking statistics.
Needs pandas, matplotlib and logomaker.
"""
import argparse, os, sys, math
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memelib import read_meme                      # noqa: E402

COLORS = {"A": "#4DAF4A", "C": "#386CB0", "G": "#FFE52C", "T": "#E41A1C"}   # sampled from the original panel
CIRCLE = "#F175AD"


def die(msg):
    sys.exit(f"[logos] ERROR: {msg}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crossrun", required=True, help="crossrun output folder")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--rank-by", default="detected", choices=["detected", "prevalence", "denovo", "sea"])
    ap.add_argument("--samples", help="comma-separated sample_ids to rank on")
    ap.add_argument("--families", help="comma-separated family_ids to draw (overrides the ranking)")
    ap.add_argument("--ncol", type=int, default=1)
    ap.add_argument("--positions", type=int, default=30,
                    help="x-axis room in motif positions (default 30, about the width of the original logos); every logo "
                         "uses max(this, widest motif) so letters are the same width in all logos")
    ap.add_argument("--width", type=float, default=7.5, help="figure width, inches")
    ap.add_argument("--row-height", type=float, default=0.9, help="height of one row of logos, inches")
    ap.add_argument("--font", type=float, default=7, help="axis text size, points")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--out", default="family_logos.png")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        import logomaker
    except ImportError:
        die("logomaker is not installed (pip install logomaker, or add it to the conda environment)")

    pf = os.path.join(a.crossrun, "tables", "table02_motif_prevalence.tsv")
    mf = os.path.join(a.crossrun, "work", "motifs", "family_reps.meme")
    for f in (pf, mf):
        if not os.path.exists(f):
            die(f"not found: {f} (run steps 02-09 first)")
    prev = pd.read_csv(pf, sep="\t")
    bg, motifs = read_meme(mf)
    reps = {m.id: m for m in motifs}

    if a.samples:
        want = [s.strip() for s in a.samples.split(",") if s.strip()]
        miss = sorted(set(want) - set(prev.sample_id))
        if miss:
            die(f"samples not in {pf}: {', '.join(miss)}")
        prev = prev[prev.sample_id.isin(want)]
    tf = lambda s: s.astype(str).str.upper().eq("TRUE")
    stats = prev.groupby("family_id").agg(
        consensus=("pwm_consensus", "first"),
        n_samples=("sample_id", "nunique"),
        samples_detected=("detected", lambda s: int(tf(s).sum())),
        samples_denovo=("discovered_denovo", lambda s: int(tf(s).sum())),
        median_pct_peaks=("pct_peaks", "median"),
        min_sea_qvalue=("sea_qvalue", "min"),
    ).reset_index()

    if a.families:
        fams = [f.strip() for f in a.families.split(",") if f.strip()]
        miss = [f for f in fams if f not in reps]
        if miss:
            die(f"families not in {mf}: {', '.join(miss)}")
        chosen = stats.set_index("family_id").reindex(fams).reset_index()
        chosen["family_id"] = fams
        rank_note = "as given"
    else:
        key = {"detected": (["samples_detected", "median_pct_peaks"], [False, False]),
               "prevalence": (["median_pct_peaks", "samples_detected"], [False, False]),
               "denovo": (["samples_denovo", "samples_detected"], [False, False]),
               "sea": (["min_sea_qvalue", "samples_detected"], [True, False])}[a.rank_by]
        stats = stats[stats.family_id.isin(reps)]
        chosen = stats.sort_values(key[0] + ["family_id"], ascending=key[1] + [True]).head(a.top)
        rank_note = f"ranked by {a.rank_by}"
        if len(chosen) < a.top:
            print(f"[logos] only {len(chosen)} families available", file=sys.stderr)
    if chosen.empty:
        die("no families to draw")
    chosen = chosen.reset_index(drop=True)
    chosen.insert(0, "panel_number", range(1, len(chosen) + 1))

    # information content (bits) of each position against the MEME file's background
    bgv = np.array([bg.get(b, 0.25) for b in "ACGT"]) if isinstance(bg, dict) else np.full(4, 0.25)
    n = len(chosen)
    span = max(a.positions, max(len(reps[f].rows) for f in chosen.family_id))   # same letter width in every logo
    nrow = math.ceil(n / a.ncol)
    fig, axes = plt.subplots(nrow, a.ncol, figsize=(a.width, a.row_height * nrow), squeeze=False)
    plt.rcParams.update({"font.size": a.font, "font.family": "sans-serif"})
    for k, ax in enumerate(axes.flat):
        if k >= n:
            ax.axis("off")
            continue
        fam = chosen.family_id[k]
        m = reps[fam]
        p = np.array(m.rows, dtype=float)
        p = p / p.sum(axis=1, keepdims=True)
        ic = (p * np.log2(np.where(p > 0, p, 1) / bgv)).sum(axis=1)          # relative entropy per position
        heights = pd.DataFrame(p * ic[:, None], columns=list("ACGT"))
        heights[heights < 0] = 0
        logomaker.Logo(heights, ax=ax, color_scheme=COLORS, font_name="DejaVu Sans", vpad=0.02, width=0.92, baseline_width=0)
        ax.set_ylim(0, 2)
        ax.set_yticks([0, 2])
        ax.set_ylabel("bits", fontsize=a.font, labelpad=1)
        ax.set_xticks([])
        ax.set_xlim(-0.5, span - 0.5)
        for sp in ("top", "right", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.plot([-0.5, len(heights) - 0.5], [0, 0], color="black", lw=0.6, clip_on=False)   # baseline under the motif only
        ax.tick_params(axis="y", labelsize=a.font, length=2, pad=1)
        # numbered pink circle to the left of the axis, as in the original panel
        ax.annotate(str(k + 1), xy=(0, 0.5), xycoords="axes fraction", xytext=(-a.font * 4.2, 0), textcoords="offset points",
                    ha="center", va="center", fontsize=a.font * 1.15, fontweight="bold", color="white",
                    bbox=dict(boxstyle="circle,pad=0.25", fc=CIRCLE, ec="none"))
    fig.subplots_adjust(left=0.09 if a.ncol == 1 else 0.07, right=0.99, top=0.97, bottom=0.03, wspace=0.35, hspace=0.35)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=a.dpi, facecolor="white")
    tsv = os.path.splitext(a.out)[0] + ".tsv"
    chosen.assign(selection=rank_note, motif_width=[len(reps[f].rows) for f in chosen.family_id]).to_csv(tsv, sep="\t", index=False)
    print(f"[logos] {n} families ({rank_note}) -> {a.out}; table -> {tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
