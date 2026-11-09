#!/usr/bin/env python3
"""Composite supplemental figures 1-3 from an existing chipseq_crossrun output folder.

  Supplemental 1  library quality (a), mapping destination (b), read-level repeat content (c)
  Supplemental 2  peak density per chromosome (a, b), satellite composition (c), enrichment (d)
  Supplemental 3  motif-family periodicity (a-e)

Everything is read from <crossrun>/tables/*.tsv and <crossrun>/work/, written by the pipeline
steps listed under --sources.  Nothing is recomputed except the re-binning of the spacing
distribution for Supplemental 3b/3c, which follows R/10_motifs.R exactly (see HARMONIC_WINDOW).

Usage:
  python3 extra/plot_supplementals.py --crossrun <outdir> [options]

Options:
  --crossrun DIR       crossrun output folder (required)
  --out DIR            output folder (default <crossrun>/figures/supplemental)
  --figures 1,2,3      which supplementals to draw (default: all)
  --qc-table FILE      TSV from extra/compile_qc_table.py (default <crossrun>/tables/supplemental2_qc_table.tsv,
                       then supp_qc_table.tsv; Supplemental 1a is skipped with a warning if absent)
  --labels-file FILE   TSV sample_id<TAB>label giving IP sample row order and figure labels
                       (default: work/samples.tsv order, sample_id used as the label)
  --colour-by COL      manifest column used to colour sample points (default virus; "" = single colour)
  --families IDS       comma-separated family ids for Supplemental 3 (default: the rule below)
  --n-families N       how many families to choose automatically (default 4)
  --min-spacings N     families need at least this many site pairs to be choosable (default 30,
                       the threshold R/10_motifs.R uses before it trusts periodicity)
  --max-classes N      most repeat/CenSat classes to show in 1c, 2c and 2d (default 10, by total signal)
  --harmonics N        multiples of the modal spacing shown in 3b (default 6)
  --unit-refs SPEC     reference lines in 3e, "bp:name,bp:name" (default "5:pentamer,171:alpha-satellite monomer")
  --format png,pdf     output formats (default png,pdf)
  --dpi N              raster resolution (default 300)
  --width-mm N         figure width (default 180)
  --seed N             seed for point jitter only (default 7; recorded in the sources file)

Default family rule for Supplemental 3: from table02_repeat_likeness, keep families with
n_spacings >= --min-spacings, take the (N-1) with the most site pairs whose flag is not
"TF-like", then add the TF-like family with the most site pairs as a contrast.  The choice is
printed and written to the sources file.  Override it with --families.

Every panel's source file and row count is written to <out>/supplemental_sources.tsv.
The script exits non-zero if a required table is missing, unreadable, or yields no rows.
"""
import argparse, csv, math, os, sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# R/10_motifs.R scores periodicity as the share of spacings with (spacing %% mode) in {0, 1, mode-1},
# i.e. a +/-1 bp window centred on each multiple of the mode.  3b re-bins with the same window so the
# per-sample bars sum to the pooled periodic_frac in table02_repeat_likeness.
HARMONIC_WINDOW = 1

MM = 1 / 25.4
INK, MUTED, GRID = "#1b1b1b", "#6e6e6e", "#e4e4e4"
OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#949494", "#5c3d99"]
SAT_COLS = ["#08306b", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#deebf7",
            "#bdbdbd", "#d9d9d9", "#e8e8e8"]

SOURCES = []
WARNINGS = []


def die(msg):
    sys.exit(f"[supplementals] ERROR: {msg}")


def warn(msg):
    WARNINGS.append(msg)
    print(f"[supplementals] WARNING: {msg}", file=sys.stderr)


def log(msg):
    print(f"[supplementals] {msg}", file=sys.stderr)


def read_tsv(path, required=True, what=""):
    """Read a TSV into a list of dicts.  Missing or empty -> die (required) or None."""
    if not os.path.exists(path):
        if required:
            die(f"{what or 'table'} not found: {path}")
        return None
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        if required:
            die(f"{what or 'table'} is empty (header only): {path}")
        return None
    return rows


def note_source(panel, path, n, detail=""):
    SOURCES.append((panel, os.path.relpath(path) if path else "", n, detail))


def num(v, default=float("nan")):
    """Parse a numeric cell.  '', NA, NaN and non-numeric text all become default."""
    if v is None:
        return default
    s = str(v).strip()
    if s == "" or s.upper() in ("NA", "NAN", "NULL", "-"):
        return default
    try:
        return float(s)
    except ValueError:
        return default


def require_rows(rows, panel, why):
    if not rows:
        die(f"{panel}: {why} left no rows. Nothing would be plotted, which is a pipeline or "
            f"filter problem rather than an empty result; fix the input instead of drawing an empty panel.")


def set_style(dpi):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 6.4,
        "axes.labelsize": 6.8, "axes.titlesize": 7.2, "axes.titleweight": "bold",
        "xtick.labelsize": 6.0, "ytick.labelsize": 6.0, "legend.fontsize": 6.0,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.0, "ytick.major.size": 2.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.dpi": dpi, "savefig.dpi": dpi,
        "axes.titlepad": 3.5, "axes.labelpad": 2.0,
    })


def panel_letter(ax, s, dx=-0.09, dy=1.0):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=INK)


def source_note(fig, text, y=0.012):
    fig.text(0.008, y, text, ha="left", va="bottom", fontsize=5.3, color=MUTED, linespacing=1.5)


def save(fig, out, stem, formats):
    paths = []
    for f in formats:
        p = os.path.join(out, f"{stem}.{f}")
        fig.savefig(p)
        paths.append(p)
    plt.close(fig)
    for p in paths:
        log(f"wrote {p}")
    return paths


# ---------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------
def load_samples(cr, labels_file, colour_by):
    """IP sample order, figure labels and the colour key, with counts logged at every step."""
    spath = os.path.join(cr, "work", "samples.tsv")
    srows = read_tsv(spath, what="work/samples.tsv")
    log(f"samples.tsv: {len(srows)} rows")
    by_id = {r["sample_id"]: r for r in srows}

    if labels_file:
        lrows = read_tsv(labels_file, what="--labels-file")
        order = [r["sample_id"] for r in lrows]
        labels = {r["sample_id"]: (r.get("label") or r["sample_id"]) for r in lrows}
        missing = [s for s in order if s not in by_id]
        if missing:
            die(f"--labels-file names samples that are not in work/samples.tsv: {', '.join(missing)}")
        extra = [r["sample_id"] for r in srows if r["sample_id"] not in labels]
        if extra:
            warn(f"{len(extra)} sample(s) in samples.tsv are not in --labels-file and are left out "
                 f"of the figures: {', '.join(extra)}")
    else:
        order = [r["sample_id"] for r in srows]
        labels = {s: s for s in order}

    groups = {}
    if colour_by:
        if colour_by not in srows[0]:
            warn(f"--colour-by {colour_by!r} is not a column of work/samples.tsv "
                 f"(columns: {', '.join(srows[0].keys())}); points are drawn in one colour")
        else:
            groups = {s: (by_id[s].get(colour_by) or "none") for s in order}
    levels = sorted({v for v in groups.values()}) if groups else []
    cmap = {lv: OI[i % len(OI)] for i, lv in enumerate(levels)}
    colours = [cmap.get(groups.get(s), OI[0]) for s in order]
    log(f"samples used: {len(order)}; colour key: {colour_by or 'none'} "
        f"({', '.join(levels) if levels else 'single colour'})")
    note_source("samples", spath, len(order), f"colour_by={colour_by or 'none'}")
    return order, labels, colours, cmap, groups


def sample_legend(cmap):
    return [Line2D([], [], marker="o", ls="", ms=2.8, color=c, label=k) for k, c in cmap.items()]


def pick_classes(rows, key, value_key, max_classes):
    """Most-represented repeat/CenSat classes, by summed value across samples."""
    tot = defaultdict(float)
    for r in rows:
        v = num(r.get(value_key), 0.0)
        if math.isfinite(v):
            tot[r[key]] += v
    keep = [k for k, _ in sorted(tot.items(), key=lambda kv: -kv[1])][:max_classes]
    dropped = len(tot) - len(keep)
    if dropped > 0:
        log(f"classes: {len(tot)} present, showing {len(keep)}, {dropped} not shown "
            f"(--max-classes); no data were filtered, only the display")
    return keep


# ---------------------------------------------------------------------------------------------
# Supplemental 1
# ---------------------------------------------------------------------------------------------
QC_METRICS = [("raw_reads", "Raw reads\n(M)", 1e-6, "{:.1f}"),
              ("read_length", "Read\nlength", 1, "{:.0f}"),
              ("pct_gc", "GC\n(%)", 1, "{:.1f}"),
              ("pct_fastqc_duplicate", "Duplicate\n(%)", 1, "{:.1f}"),
              ("pct_reads_with_adapter", "Adapter\n(%)", 1, "{:.1f}"),
              ("pct_reads_passing_trimming", "Passing\ntrim (%)", 1, "{:.1f}"),
              ("pct_properly_paired", "Properly\npaired (%)", 1, "{:.1f}")]


def supp1(cr, out, args, order, labels, colours, cmap):
    tabs = os.path.join(cr, "tables")

    # ---- a: QC matrix -----------------------------------------------------------------------
    qc_path = args.qc_table
    if not qc_path:
        for cand in ("supplemental2_qc_table.tsv", "supp_qc_table.tsv"):
            if os.path.exists(os.path.join(tabs, cand)):
                qc_path = os.path.join(tabs, cand)
                break
    qc = read_tsv(qc_path, required=False, what="QC table") if qc_path else None
    if qc is None:
        warn("no QC table found (extra/compile_qc_table.py writes it); Supplemental 1a is left blank. "
             "Pass --qc-table to fill it.")
        qc_rows, qc_labels, qc_is_ip = [], [], []
    else:
        ip = [r for r in qc if (r.get("library_type") or "IP").upper() == "IP"]
        ct = [r for r in qc if (r.get("library_type") or "IP").upper() != "IP"]
        rank = {s: i for i, s in enumerate(order)}
        ip.sort(key=lambda r: rank.get(r["sample_id"], len(rank)))
        qc_rows = ip + ct
        qc_labels = [r.get("label") or labels.get(r["sample_id"], r["sample_id"]) for r in qc_rows]
        qc_is_ip = [True] * len(ip) + [False] * len(ct)
        log(f"QC table: {len(qc)} libraries ({len(ip)} IP, {len(ct)} control)")
        note_source("S1a", qc_path, len(qc_rows))

    # ---- b: contig reads --------------------------------------------------------------------
    cpath = os.path.join(tabs, "table00c_contig_reads.tsv")
    cr_rows = read_tsv(cpath, what="table00c_contig_reads")
    seen = {}
    for r in cr_rows:
        if r["sample_id"] not in order:
            continue
        lib = (r.get("library") or "IP").upper()
        if lib != "IP":
            continue
        seen.setdefault(r["sample_id"], {})
        g, contig = r.get("group", ""), r.get("contig", "")
        pct = num(r.get("pct_reads"), 0.0)
        if not math.isfinite(pct):
            pct = 0.0
        key = "Host nuclear" if g == "nuclear" else ("Unaligned" if g == "unaligned" else contig)
        seen[r["sample_id"]][key] = seen[r["sample_id"]].get(key, 0.0) + pct
    require_rows(seen, "S1b", "keeping IP libraries of the listed samples")
    missing_b = [s for s in order if s not in seen]
    if missing_b:
        warn(f"S1b: no table00c rows for {len(missing_b)} sample(s): {', '.join(missing_b)}; "
             f"they are drawn empty rather than dropped")
    keys = set()
    for d in seen.values():
        keys |= set(d)
    small = sorted(k for k in keys if k not in ("Host nuclear", "Unaligned"))
    stack_keys = ["Host nuclear"] + small + ["Unaligned"]
    log(f"S1b: {len(seen)} IP libraries, contig groups: {', '.join(stack_keys)}")
    note_source("S1b", cpath, len(cr_rows), f"stack={'|'.join(stack_keys)}")

    # ---- c: reads in repeats ----------------------------------------------------------------
    rpath = os.path.join(tabs, "table03d_reads_in_repeats.tsv")
    rr = read_tsv(rpath, required=False, what="table03d_reads_in_repeats")
    if rr is None:
        warn("table03d_reads_in_repeats.tsv not found (set repeats.reads_in_repeats: true and rerun "
             "step 07); Supplemental 1c is left blank.")
        rr, classes = [], []
    else:
        rr = [r for r in rr if r["sample_id"] in order]
        require_rows(rr, "S1c", "keeping the listed samples")
        for r in rr:
            r["_label"] = r.get("label_name") or r.get("label") or "?"
        classes = pick_classes(rr, "_label", "ip_frac", args.max_classes)
        log(f"S1c: {len(rr)} sample x class rows, {len(classes)} classes shown")
        note_source("S1c", rpath, len(rr), f"classes={'|'.join(classes)}")

    # ---- draw --------------------------------------------------------------------------------
    W = args.width_mm * MM
    n = max(len(qc_rows), len(order))
    fig = plt.figure(figsize=(W, 205 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.62], width_ratios=[1.35, 1.0],
                          left=0.205, right=0.975, top=0.950, bottom=0.120, hspace=0.42, wspace=0.17)

    ax = fig.add_subplot(gs[0, 0])
    if qc_rows:
        vals = np.array([[num(r.get(k)) * sc for k, _, sc, _ in QC_METRICS] for r in qc_rows])
        mu, sd = np.nanmean(vals, 0), np.nanstd(vals, 0)
        z = (vals - mu) / np.where(sd > 0, sd, 1)
        ax.imshow(np.nan_to_num(z), cmap="RdYlBu_r", vmin=-2.2, vmax=2.2, aspect="auto")
        for i in range(len(qc_rows)):
            for j, (_, _, _, fmt) in enumerate(QC_METRICS):
                v = vals[i, j]
                ax.text(j, i, "" if not math.isfinite(v) else fmt.format(v),
                        ha="center", va="center", fontsize=5.1, color=INK)
        ax.set_xticks(range(len(QC_METRICS)))
        ax.set_xticklabels([m[1] for m in QC_METRICS], fontsize=5.6)
        ax.set_yticks(range(len(qc_rows))); ax.set_yticklabels(qc_labels)
        for t, isip, s in zip(ax.get_yticklabels(), qc_is_ip, [r["sample_id"] for r in qc_rows]):
            t.set_color(dict(zip(order, colours)).get(s, MUTED) if isip else MUTED)
        ax.set_xticks(np.arange(-.5, len(QC_METRICS), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(qc_rows), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.8)
        ax.tick_params(which="minor", length=0); ax.tick_params(axis="both", length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        if any(qc_is_ip) and not all(qc_is_ip):
            ax.axhline(sum(qc_is_ip) - 0.5, color=INK, lw=0.9)
        axk = ax.inset_axes([0.0, -0.165, 0.30, 0.020])
        axk.imshow(np.linspace(-2.2, 2.2, 128)[None, :], cmap="RdYlBu_r", aspect="auto")
        axk.set_xticks([0, 127]); axk.set_xticklabels(["low", "high"], fontsize=5.0)
        axk.set_yticks([]); axk.tick_params(length=0, pad=1)
        for sp in axk.spines.values():
            sp.set_linewidth(0.4)
        axk.set_title("colour: z-score within column", fontsize=5.2, color=MUTED,
                      fontweight="normal", pad=2, loc="left")
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no QC table\n(--qc-table)", ha="center", va="center",
                fontsize=6.5, color=MUTED, transform=ax.transAxes)
    ax.set_title("Library quality summary", loc="left")
    panel_letter(ax, "a", dx=-0.30)

    gsb = gs[0, 1].subgridspec(1, 2, width_ratios=[1.0, 0.42], wspace=0.10)
    ax = fig.add_subplot(gsb[0, 0])
    yb = np.arange(len(order))
    left = np.zeros(len(order))
    pal = {"Host nuclear": "#b9c6d4", "Unaligned": "#d9d9d9"}
    for i, k in enumerate(small):
        pal[k] = OI[(i + 3) % len(OI)]
    for k in stack_keys:
        v = np.array([seen.get(s, {}).get(k, 0.0) for s in order])
        ax.barh(yb, v, left=left, color=pal[k], height=0.74, label=k, edgecolor="white", linewidth=0.4)
        left += v
    ax.set_xlim(0, 100); ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 25, 50, 75, 100]); ax.set_xlabel("Reads (% of library)")
    ax.set_title("Mapping destination", loc="left")
    ax.legend(ncol=2, loc="upper left", bbox_to_anchor=(-0.02, -0.105), handlelength=0.8,
              handleheight=0.8, columnspacing=0.7, fontsize=5.6, borderpad=0)
    panel_letter(ax, "b", dx=-0.135)

    ax2 = fig.add_subplot(gsb[0, 1])
    if small:
        h = 0.74 / len(small)
        for i, k in enumerate(small):
            v = np.array([max(seen.get(s, {}).get(k, 0.0), 1.5e-3) for s in order])
            ax2.barh(yb + (i - (len(small) - 1) / 2) * h, v, color=pal[k], height=h * 0.9)
        ax2.set_xscale("log"); ax2.set_xlim(1.5e-3, 20)
        ax2.set_xticks([1e-2, 1e0])
    ax2.set_ylim(len(order) - 0.5, -0.5)
    ax2.set_yticks([]); ax2.spines["left"].set_visible(False)
    ax2.set_xlabel("% of library\n(log)")
    ax2.set_axisbelow(True); ax2.xaxis.grid(True, color=GRID, lw=0.5)
    ax2.set_title("chrM / viral\ndetail", loc="left", fontsize=5.8, color=MUTED, fontweight="normal")

    ax = fig.add_subplot(gs[1, :])
    if classes:
        x = np.arange(len(classes))
        byc = defaultdict(lambda: ([], []))
        for r in rr:
            if r["_label"] in classes:
                byc[r["_label"]][0].append(num(r.get("ip_frac")) * 100)
                byc[r["_label"]][1].append(num(r.get("control_frac")) * 100)
        wd = 0.34
        ipm = [np.nanmean(byc[c][0]) if byc[c][0] else np.nan for c in classes]
        ctm = [np.nanmean(byc[c][1]) if byc[c][1] else np.nan for c in classes]
        ax.bar(x - wd / 2, ipm, wd, color=OI[0], alpha=0.30, edgecolor=OI[0], linewidth=0.6, label="IP (mean)")
        ax.bar(x + wd / 2, ctm, wd, color=OI[6], alpha=0.30, edgecolor=OI[6], linewidth=0.6, label="Control (mean)")
        rng = np.random.default_rng(args.seed)
        ccol = dict(zip(order, colours))
        for j, c in enumerate(classes):
            ipv = [v for v in byc[c][0] if math.isfinite(v)]
            ids = [r["sample_id"] for r in rr if r["_label"] == c and math.isfinite(num(r.get("ip_frac")))]
            ax.scatter(np.full(len(ipv), x[j] - wd / 2) + rng.uniform(-0.09, 0.09, len(ipv)), ipv,
                       s=3.6, c=[ccol.get(s, OI[0]) for s in ids], linewidths=0, zorder=3)
            ctv = [v for v in byc[c][1] if math.isfinite(v)]
            ax.scatter(np.full(len(ctv), x[j] + wd / 2) + rng.uniform(-0.06, 0.06, len(ctv)), ctv,
                       s=3.6, c=INK, linewidths=0, zorder=3)
        ax.set_yscale("log")
        ax.set_xticks(x); ax.set_xticklabels(classes, rotation=28, ha="right")
        ax.set_ylabel("Mapped reads in class (%, log)")
        ax.set_axisbelow(True); ax.yaxis.grid(True, color=GRID, lw=0.5)
        ax.legend(handles=sample_legend(cmap) + [Line2D([], [], marker="o", ls="", ms=2.8, color=INK, label="control")]
                  + ax.get_legend_handles_labels()[0], ncol=5, loc="upper right", handlelength=0.9)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "table03d_reads_in_repeats.tsv not found\n"
                          "(set repeats.reads_in_repeats: true and rerun step 07)",
                ha="center", va="center", fontsize=6.5, color=MUTED, transform=ax.transAxes)
    ax.set_title("Read-level repeat and satellite content, IP versus matched control", loc="left")
    panel_letter(ax, "c", dx=-0.095)

    source_note(fig,
                "a  extra/compile_qc_table.py output (raw_reads, read_length, pct_gc, pct_fastqc_duplicate, "
                "pct_reads_with_adapter, pct_reads_passing_trimming, pct_properly_paired)\n"
                "b  tables/table00c_contig_reads.tsv — pct_reads by contig group (nuclear / excluded / viral / unaligned)\n"
                "c  tables/table03d_reads_in_repeats.tsv — ip_frac and control_frac per repeat and CenSat class (step 07)")
    return save(fig, out, "supplemental1_library_qc", args.format)


# ---------------------------------------------------------------------------------------------
# Supplemental 2
# ---------------------------------------------------------------------------------------------
def supp2(cr, out, args, order, labels, colours, cmap):
    tabs = os.path.join(cr, "tables")
    ppath = os.path.join(tabs, "table00b_peaks_by_chromosome.tsv")
    pc = read_tsv(ppath, what="table00b_peaks_by_chromosome")
    pc = [r for r in pc if r["sample_id"] in order]
    require_rows(pc, "S2a", "keeping the listed samples")
    chroms = []
    for r in pc:
        if r["chrom"] not in chroms:
            chroms.append(r["chrom"])
    clen = {r["chrom"]: num(r.get("chrom_length")) for r in pc}
    dens = np.full((len(order), len(chroms)), np.nan)
    npk = np.zeros(len(chroms))
    ri = {s: i for i, s in enumerate(order)}
    ci = {c: i for i, c in enumerate(chroms)}
    for r in pc:
        dens[ri[r["sample_id"]], ci[r["chrom"]]] = num(r.get("peaks_per_mb"))
        npk[ci[r["chrom"]]] += num(r.get("n_peaks"), 0.0)
    log(f"S2a: {len(pc)} sample x chromosome rows, {len(order)} samples, {len(chroms)} chromosomes, "
        f"{int(npk.sum())} peaks pooled")
    note_source("S2a/b", ppath, len(pc))

    rpath = os.path.join(tabs, "table03_repeat_enrichment.tsv")
    re_rows = read_tsv(rpath, what="table03_repeat_enrichment")
    before = len(re_rows)
    re_rows = [r for r in re_rows if r["sample_id"] in order and r.get("feature_set") == "censat"]
    log(f"S2c/d: table03_repeat_enrichment {before} rows -> {len(re_rows)} after keeping "
        f"feature_set=censat and the listed samples")
    require_rows(re_rows, "S2c/d", "keeping feature_set=censat and the listed samples")
    for r in re_rows:
        r["_label"] = r.get("label_name") or r.get("label") or "?"
    classes = pick_classes(re_rows, "_label", "pct_of_peaks", args.max_classes)
    note_source("S2c/d", rpath, len(re_rows), f"feature_set=censat; classes={'|'.join(classes)}")

    W = args.width_mm * MM
    n = len(order)
    fig = plt.figure(figsize=(W, 205 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.72, 0.80], width_ratios=[1, 1],
                          left=0.185, right=0.965, top=0.950, bottom=0.115, hspace=0.95, wspace=0.42)

    ax = fig.add_subplot(gs[0, :])
    im = ax.imshow(np.nan_to_num(dens), cmap="magma_r", aspect="auto", vmin=0,
                   vmax=np.nanmax(dens) if np.isfinite(dens).any() else 1)
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=5.4)
    ax.set_yticks(range(n)); ax.set_yticklabels([labels[s] for s in order])
    for t, c in zip(ax.get_yticklabels(), colours):
        t.set_color(c)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("Chromosome")
    ax.set_title("Peak density per chromosome", loc="left")
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.022)
    cb.set_label("Peaks per Mb", fontsize=5.8); cb.ax.tick_params(labelsize=5.2, length=1.5)
    cb.outline.set_linewidth(0.4)
    panel_letter(ax, "a", dx=-0.24)

    ax = fig.add_subplot(gs[1, 0])
    L = np.array([clen.get(c, np.nan) for c in chroms])
    exp = npk.sum() * L / np.nansum(L)
    x = np.arange(len(chroms))
    ax.bar(x, npk, 0.74, color=OI[0], alpha=0.75, label="Observed", linewidth=0)
    ax.plot(x, exp, color=OI[1], lw=0.9, ls="--", label="Expected from length")
    ax.set_xticks(x); ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=5.0)
    ax.set_xlabel("Chromosome"); ax.set_ylabel("Peaks (all samples)")
    ax.set_axisbelow(True); ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.legend(loc="upper right", handlelength=1.1)
    ax.set_title("Pooled peaks vs chromosome length", loc="left")
    panel_letter(ax, "b", dx=-0.22)

    ax = fig.add_subplot(gs[1, 1])
    pct = {(r["sample_id"], r["_label"]): num(r.get("pct_of_peaks"), 0.0) for r in re_rows}
    yb = np.arange(n); left = np.zeros(n)
    for i, c in enumerate(classes):
        v = np.array([pct.get((s, c), 0.0) for s in order])
        ax.barh(yb, v, left=left, color=SAT_COLS[i % len(SAT_COLS)], height=0.74, label=c,
                edgecolor="white", linewidth=0.35)
        left += v
    rest = 100 - left
    if np.any(rest > 0.5):
        ax.barh(yb, np.clip(rest, 0, None), left=left, color="#e8e8e8", height=0.74,
                label="other / none", edgecolor="white", linewidth=0.35)
    ax.set_xlim(0, 100); ax.set_ylim(n - 0.5, -0.5)
    ax.set_yticks(yb); ax.set_yticklabels([labels[s] for s in order], fontsize=5.0)
    for t, c in zip(ax.get_yticklabels(), colours):
        t.set_color(c)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Peaks in class (%)")
    ax.set_title("Satellite composition of peaks", loc="left")
    ax.legend(ncol=2, loc="upper left", bbox_to_anchor=(0, -0.20), handlelength=0.9,
              handleheight=0.9, columnspacing=0.8, fontsize=5.6, borderpad=0)
    panel_letter(ax, "c", dx=-0.30)

    ax = fig.add_subplot(gs[2, :])
    ccol = dict(zip(order, colours))
    rng = np.random.default_rng(args.seed)
    for j, c in enumerate(classes):
        vs = [(r["sample_id"], num(r.get("log2_oe"))) for r in re_rows if r["_label"] == c]
        vs = [(s, v) for s, v in vs if math.isfinite(v)]
        if not vs:
            continue
        xs = np.array([v for _, v in vs])
        ax.scatter(xs, np.full(len(xs), j) + rng.uniform(-0.16, 0.16, len(xs)), s=7,
                   c=[ccol.get(s, OI[0]) for s, _ in vs], linewidths=0, zorder=3)
        lo, hi = np.percentile(xs, [2.5, 97.5]) if len(xs) > 1 else (xs[0], xs[0])
        ax.hlines(j, lo, hi, color=INK, lw=1.4, zorder=2)
        ax.plot(np.median(xs), j, "|", color=INK, ms=7, mew=1.4, zorder=4)
    ax.axvline(0, color=MUTED, lw=0.7, ls=":")
    ax.set_yticks(range(len(classes))); ax.set_yticklabels(classes)
    ax.set_ylim(len(classes) - 0.4, -0.6)
    ax.set_xlabel("log$_2$ (observed / expected peaks)")
    ax.set_axisbelow(True); ax.xaxis.grid(True, color=GRID, lw=0.5)
    ax.set_title("Peak enrichment by satellite class (bar = sample range, tick = median)", loc="left")
    ax.legend(handles=sample_legend(cmap), ncol=max(1, len(cmap)), loc="lower right", handlelength=0.9)
    panel_letter(ax, "d", dx=-0.21)

    source_note(fig,
                "a, b  tables/table00b_peaks_by_chromosome.tsv — n_peaks, chrom_length, peaks_per_mb\n"
                "c  tables/table03_repeat_enrichment.tsv (feature_set = censat) — pct_of_peaks\n"
                "d  tables/table03_repeat_enrichment.tsv — log2_oe against the configured null "
                "(repeats.null_model, repeats.n_shuffles); see null_* columns in that table")
    return save(fig, out, "supplemental2_peak_distribution", args.format)


# ---------------------------------------------------------------------------------------------
# Supplemental 3
# ---------------------------------------------------------------------------------------------
def load_spacing(cr, order):
    """work/scan/<sample>.peaks.spacing.tsv -> pooled[fam][spacing] and per_sample[fam][sample][spacing]."""
    sdir = os.path.join(cr, "work", "scan")
    if not os.path.isdir(sdir):
        die(f"work/scan not found: {sdir} (run step 05 first)")
    files = sorted(f for f in os.listdir(sdir) if f.endswith(".peaks.spacing.tsv"))
    if not files:
        die(f"no *.peaks.spacing.tsv in {sdir} (step 05 writes them for the peaks set)")
    pooled = defaultdict(lambda: defaultdict(float))
    per = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    kept = skipped = 0
    for f in files:
        for r in read_tsv(os.path.join(sdir, f), what=f):
            s = r["sample_id"]
            if s not in order:
                skipped += 1
                continue
            d, nn = int(num(r["spacing"], 0)), num(r["n"], 0.0)
            if d < 2:          # R/10_motifs.R drops spacing < 2 before taking the mode
                continue
            pooled[r["family_id"]][d] += nn
            per[r["family_id"]][s][d] += nn
            kept += 1
    log(f"S3: {len(files)} spacing files, {kept} rows kept, {skipped} rows from samples not listed")
    note_source("S3a/b/c", sdir, kept, f"{len(files)} files; spacing >= 2 as in R/10_motifs.R")
    return pooled, per


def mode_of(hist):
    return max(hist.items(), key=lambda kv: kv[1])[0] if hist else None


def harmonic_shares(hist, mode, nh):
    """Share of pairs within +/-HARMONIC_WINDOW bp of k x mode, k = 1..nh (the periodic_frac window)."""
    tot = sum(hist.values())
    if not tot or not mode:
        return np.full(nh, np.nan)
    return np.array([sum(v for d, v in hist.items() if abs(d - k * mode) <= HARMONIC_WINDOW) / tot
                     for k in range(1, nh + 1)])


def choose_families(rl, args):
    rows = [r for r in rl if num(r.get("n_spacings"), 0) >= args.min_spacings]
    if not rows:
        die(f"no motif family has at least --min-spacings {args.min_spacings} site pairs in "
            f"table02_repeat_likeness; pass --families explicitly or lower the threshold")
    rows.sort(key=lambda r: -num(r.get("n_spacings"), 0))
    per = [r for r in rows if (r.get("flag") or "") != "TF-like"]
    tf = [r for r in rows if (r.get("flag") or "") == "TF-like"]
    chosen = per[:max(args.n_families - 1, 0)] + tf[:1]
    if len(chosen) < args.n_families:
        for r in rows:
            if r not in chosen and len(chosen) < args.n_families:
                chosen.append(r)
    chosen = chosen[:args.n_families]
    log("S3 families chosen automatically (rule: most site pairs, not TF-like, plus the top TF-like "
        "family as a contrast). Override with --families:")
    for r in chosen:
        log(f"    {r['family_id']}  flag={r.get('flag')}  n_spacings={r.get('n_spacings')}  "
            f"spacing_mode={r.get('spacing_mode')}")
    return [r["family_id"] for r in chosen]


def supp3(cr, out, args, order, labels, colours, cmap):
    tabs = os.path.join(cr, "tables")
    rl_path = os.path.join(tabs, "table02_repeat_likeness.tsv")
    rl = read_tsv(rl_path, what="table02_repeat_likeness")
    byfam = {r["family_id"]: r for r in rl}
    pooled, per = load_spacing(cr, order)

    if args.families:
        fams = [f.strip() for f in args.families.split(",") if f.strip()]
        missing = [f for f in fams if f not in byfam]
        if missing:
            die(f"--families not in table02_repeat_likeness: {', '.join(missing)}")
    else:
        fams = choose_families(rl, args)
    note_source("S3", rl_path, len(rl), f"families={'|'.join(fams)}"
                + ("" if args.families else f"; auto, min_spacings={args.min_spacings}"))

    pv_path = os.path.join(tabs, "table02_motif_prevalence.tsv")
    pv = read_tsv(pv_path, what="table02_motif_prevalence")
    pv = [r for r in pv if r["family_id"] in fams and r["sample_id"] in order]
    require_rows(pv, "S3d", "keeping the chosen families and listed samples")
    note_source("S3d", pv_path, len(pv))

    cols = [OI[i % len(OI)] for i in range(len(fams))]
    flabs = []
    for f in fams:
        top = byfam[f].get("top_repeat_label") or byfam[f].get("flag") or ""
        flabs.append(f"{f}  {top}" if top else f)

    W = args.width_mm * MM
    nIP = len(order)
    fig = plt.figure(figsize=(W, 240 * MM))
    gsA = fig.add_gridspec(2, len(fams), left=0.105, right=0.985, top=0.955, bottom=0.600,
                           hspace=0.82, wspace=0.42)
    gsB = fig.add_gridspec(1, 2, left=0.215, right=0.985, top=0.505, bottom=0.355, wspace=0.62)
    y = np.arange(len(fams))
    ccol = dict(zip(order, colours))
    rng = np.random.default_rng(args.seed)

    def famy(ax):
        ax.set_yticks(y); ax.set_yticklabels(flabs, fontsize=5.8)
        ax.set_ylim(len(fams) - 0.45, -0.55)
        ax.tick_params(axis="y", length=0)
        ax.set_axisbelow(True); ax.xaxis.grid(True, color=GRID, lw=0.5)

    # The pooled mode and the periodicity verdict come from table02_repeat_likeness, so the figure
    # always agrees with the table. R/10_motifs.R leaves periodic_frac NA when it will not trust the
    # mode (mode < 3), and a decaying histogram always has a "mode" at its smallest spacing, so a
    # family with NA periodic_frac is drawn as aperiodic rather than given a meaningless period.
    modes, periodic = {}, {}
    for f in fams:
        tab = num(byfam[f].get("spacing_mode"))
        pf = num(byfam[f].get("periodic_frac"))
        periodic[f] = math.isfinite(pf)
        if math.isfinite(tab) and tab > 0:
            modes[f] = int(tab)
        else:
            modes[f] = mode_of(pooled.get(f, {}))
            if periodic[f]:
                warn(f"S3: no spacing_mode in table02_repeat_likeness for {f}; using the mode "
                     f"recomputed from work/scan ({modes[f]})")
        if not periodic[f]:
            modes[f] = None
    log("S3 pooled modes (from table02_repeat_likeness.spacing_mode): "
        + ", ".join(f"{f}={modes[f] if modes[f] else 'aperiodic'}" for f in fams))

    # a: spacing histogram
    for k, f in enumerate(fams):
        ax = fig.add_subplot(gsA[0, k])
        hist = pooled.get(f, {})
        if not hist:
            ax.axis("off")
            ax.text(0.5, 0.5, f"{f}\nno spacings", ha="center", va="center", fontsize=6, color=MUTED,
                    transform=ax.transAxes)
            continue
        mode = modes[f]
        xmax = min(max(hist), 220 if (mode is None or mode > 40) else 70)
        d = np.arange(2, xmax + 1)
        h = np.array([hist.get(int(i), 0.0) for i in d])
        ax.fill_between(d, h, color=cols[k], alpha=0.85, linewidth=0)
        if mode:
            for mm in range(1, min(int(xmax // mode), 8) + 1):
                ax.axvline(mode * mm, color=INK, lw=0.45, ls=(0, (2, 2.5)), alpha=0.6)
            ax.text(0.97, 0.92, f"mode {mode} bp", transform=ax.transAxes, fontsize=5.6,
                    ha="right", va="top", color=INK)
        else:
            ax.text(0.97, 0.92, "no periodic mode", transform=ax.transAxes, fontsize=5.6,
                    ha="right", va="top", color=MUTED, style="italic")
        ax.set_xlim(2, xmax); ax.set_ylim(0, (h.max() if h.size else 1) * 1.20)
        ax.set_xlabel("Spacing (bp)")
        if k == 0:
            ax.set_ylabel("Adjacent site pairs"); panel_letter(ax, "a", dx=-0.46)
        ax.set_title(f"{f}\n{byfam[f].get('flag') or ''}", loc="left", fontsize=6.2)
        ax.set_axisbelow(True); ax.yaxis.grid(True, color=GRID, lw=0.5)

    # b: harmonic structure
    ks = np.arange(1, args.harmonics + 1)
    for k, f in enumerate(fams):
        ax = fig.add_subplot(gsA[1, k])
        mode = modes[f]
        shares = harmonic_shares(pooled.get(f, {}), mode, args.harmonics)
        ax.bar(ks, np.nan_to_num(shares), 0.66, color=cols[k], alpha=0.32, edgecolor=cols[k],
               linewidth=0.7, zorder=2)
        for s in order:
            hs = harmonic_shares(per.get(f, {}).get(s, {}), mode, args.harmonics)
            good = np.isfinite(hs)
            ax.scatter(ks[good] + rng.uniform(-0.17, 0.17, good.sum()), hs[good], s=3.0,
                       c=ccol.get(s, OI[0]), linewidths=0, zorder=4)
        ax.set_xticks(ks); ax.set_xticklabels([f"{i}×" for i in ks])
        ax.set_xlim(0.45, args.harmonics + 0.55)
        top = np.nanmax(shares) if np.isfinite(shares).any() else 1
        ax.set_ylim(0, top * 1.45)
        ax.set_xlabel("Multiple of modal spacing")
        pf = num(byfam[f].get("periodic_frac"))
        ax.text(0.97, 0.93, f"periodic_frac {pf:.2f}" if math.isfinite(pf) else "periodic_frac n/a",
                transform=ax.transAxes, fontsize=5.6, ha="right", va="top",
                color=INK if math.isfinite(pf) else MUTED)
        if k == 0:
            ax.set_ylabel("Share of site pairs"); panel_letter(ax, "b", dx=-0.46)
            ax.legend(handles=sample_legend(cmap), loc="upper left", bbox_to_anchor=(-0.02, -0.32),
                      ncol=max(1, len(cmap)), handlelength=1.3, fontsize=5.6, columnspacing=1.2, borderpad=0)
        ax.set_axisbelow(True); ax.yaxis.grid(True, color=GRID, lw=0.5)

    # c: per-sample mode
    ax = fig.add_subplot(gsB[0, 0])
    for i, f in enumerate(fams):
        ms, ss = [], []
        for s in order:
            m = mode_of(per.get(f, {}).get(s, {}))
            if m:
                ms.append(m); ss.append(s)
        if ms:
            ax.scatter(ms, np.full(len(ms), i) + rng.uniform(-0.17, 0.17, len(ms)), s=7,
                       c=[ccol.get(s, OI[0]) for s in ss], linewidths=0, zorder=3)
        pooled_mode = modes[f]
        if pooled_mode and periodic[f]:
            ax.plot(pooled_mode, i, "|", color=INK, ms=9, mew=1.3, zorder=4)
            within = sum(1 for m in ms if abs(m - pooled_mode) <= 2)
            note = f"{within}/{len(ms)} within ±2 bp"
        else:
            note = "no mode"
        ax.text(0.985, i - 0.42, note, transform=ax.get_yaxis_transform(), fontsize=5.2,
                color=INK, ha="right", va="center")
    ax.set_xscale("log"); famy(ax)
    ax.set_xlabel("Modal spacing per sample (bp, log)")
    ax.set_title("Same period, recovered in each library", loc="left", fontsize=7.0)
    panel_letter(ax, "c", dx=-0.56)

    # d: sites per peak
    ax = fig.add_subplot(gsB[0, 1])
    for i, f in enumerate(fams):
        vs = [(r["sample_id"], num(r.get("mean_sites_per_hit_peak"))) for r in pv if r["family_id"] == f]
        vs = [(s, v) for s, v in vs if math.isfinite(v) and v > 0]
        if vs:
            ax.scatter([v for _, v in vs], np.full(len(vs), i) + rng.uniform(-0.17, 0.17, len(vs)),
                       s=7, c=[ccol.get(s, OI[0]) for s, _ in vs], linewidths=0, zorder=3)
        mu = num(byfam[f].get("mean_sites_per_hit_peak"))
        if math.isfinite(mu) and mu > 0:
            ax.plot(mu, i, "|", color=INK, ms=9, mew=1.3, zorder=4)
        mx = num(byfam[f].get("max_sites_in_one_peak"))
        if not math.isfinite(mx):
            mx = max((num(r.get("max_sites")) for r in pv if r["family_id"] == f), default=float("nan"))
        if math.isfinite(mx) and mx > 0:
            ax.plot(mx, i, marker=5, ms=4.0, color=INK, ls="", zorder=4)
    if args.sites_threshold is not None:
        ax.axvline(args.sites_threshold, color=INK, lw=0.7, ls=(0, (3, 2)), alpha=0.7)
        ax.text(args.sites_threshold * 1.15, len(fams) - 0.62,
                f"repeat_min_sites_per_peak = {args.sites_threshold:g}", fontsize=5.2, color=INK,
                va="center", ha="left")
    ax.set_xscale("log"); famy(ax)
    ax.set_xlabel("Sites per peak that has a site (log)")
    ax.set_title("Multiple copies within one peak", loc="left", fontsize=7.0)
    ax.legend(handles=[Line2D([], [], marker="|", ls="", ms=7, mew=1.3, color=INK, label="pooled mean"),
                       Line2D([], [], marker=5, ls="", ms=4, color=INK, label="max in one peak")],
              loc="upper left", bbox_to_anchor=(0.0, -0.33), ncol=2, handlelength=1.1,
              fontsize=5.6, columnspacing=1.2, borderpad=0)
    panel_letter(ax, "d", dx=-0.42)

    # e: predicted gap vs annotated units
    ax = fig.add_axes([0.245, 0.135, 0.515, 0.125])
    for i, f in enumerate(fams):
        m = modes[f]
        if not m or not periodic[f]:
            ax.text(0.17, i, "no periodic spacing detected", transform=ax.get_yaxis_transform(),
                    fontsize=5.6, color=MUTED, va="center", style="italic")
            continue
        spread = [mode_of(per.get(f, {}).get(s, {})) for s in order]
        spread = [v for v in spread if v]
        err = (np.std(spread) if len(spread) > 1 else 0.0)
        ax.errorbar(m, i, xerr=err, fmt="o", ms=4.2, color=cols[i], ecolor=cols[i],
                    elinewidth=1.1, capsize=1.8, zorder=3)
    for ref, txt in args.unit_refs:
        ax.axvline(ref, color=INK, lw=0.6, ls=(0, (3, 2)), alpha=0.55)
        ax.text(ref * 1.08, -1.05, txt, fontsize=5.3, color=INK, va="bottom", ha="left")
    ax.set_xscale("log")
    ax.set_yticks(y); ax.set_yticklabels(flabs, fontsize=5.8)
    ax.set_ylim(len(fams) - 0.4, -1.45)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Modal spacing (bp, log scale); bar = SD across samples")
    ax.set_axisbelow(True); ax.xaxis.grid(True, color=GRID, lw=0.5)
    ax.set_title("Repeat gap against annotated satellite unit lengths", loc="left", fontsize=7.0)
    panel_letter(ax, "e", dx=-0.30)
    for i, f in enumerate(fams):
        r = byfam[f]
        pin = num(r.get("pct_sites_in_repeats"))
        ax.text(1.03, i, f"{r.get('flag') or ''}\nn = {num(r.get('n_spacings'), 0):,.0f} pairs"
                         + (f" · {pin:.0f}% of sites in repeats" if math.isfinite(pin) else ""),
                transform=ax.get_yaxis_transform(), fontsize=5.3, color=MUTED,
                va="center", ha="left", linespacing=1.4)

    source_note(fig,
                "a  work/scan/<sample>.peaks.spacing.tsv (step 05; same-strand adjacent FIMO sites, capped at motifs.spacing_max), pooled over samples\n"
                f"b  same files, re-binned within ±{HARMONIC_WINDOW} bp of each multiple of the modal spacing — the window R/10_motifs.R uses for periodic_frac\n"
                "c  same files, one modal spacing per sample; tick = pooled spacing_mode from tables/table02_repeat_likeness\n"
                "d  tables/table02_motif_prevalence.tsv (mean_sites_per_hit_peak, max_sites per sample); pooled values from\n"
                "     tables/table02_repeat_likeness; dashed line = motifs.repeat_min_sites_per_peak\n"
                "e  pooled modal spacing with the across-sample SD; flags and counts from tables/table02_repeat_likeness",
                y=0.012)
    return save(fig, out, "supplemental3_motif_periodicity", args.format)


# ---------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crossrun", required=True)
    ap.add_argument("--out")
    ap.add_argument("--figures", default="1,2,3")
    ap.add_argument("--qc-table")
    ap.add_argument("--labels-file")
    ap.add_argument("--colour-by", default="virus")
    ap.add_argument("--families")
    ap.add_argument("--n-families", type=int, default=4)
    ap.add_argument("--min-spacings", type=int, default=30)
    ap.add_argument("--max-classes", type=int, default=10)
    ap.add_argument("--harmonics", type=int, default=6)
    ap.add_argument("--unit-refs", default="5:pentamer,171:alpha-satellite monomer")
    ap.add_argument("--sites-threshold", type=float, default=3.0,
                    help="dashed line in S3d; set to a negative value to omit (default 3, "
                         "matching motifs.repeat_min_sites_per_peak)")
    ap.add_argument("--format", default="png,pdf")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--width-mm", type=float, default=180)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    cr = os.path.abspath(a.crossrun)
    if not os.path.isdir(cr):
        die(f"--crossrun is not a directory: {cr}")
    a.out = os.path.abspath(a.out) if a.out else os.path.join(cr, "figures", "supplemental")
    os.makedirs(a.out, exist_ok=True)
    a.format = [f.strip() for f in a.format.split(",") if f.strip()]
    if a.sites_threshold is not None and a.sites_threshold < 0:
        a.sites_threshold = None
    refs = []
    for part in a.unit_refs.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            die(f"--unit-refs entry must be bp:name, got {part!r}")
        bp, name = part.split(":", 1)
        refs.append((float(bp), name.strip()))
    a.unit_refs = refs

    set_style(a.dpi)
    order, labels, colours, cmap, _ = load_samples(cr, a.labels_file, a.colour_by)
    want = {s.strip() for s in a.figures.split(",") if s.strip()}

    made = []
    if "1" in want:
        made += supp1(cr, a.out, a, order, labels, colours, cmap)
    if "2" in want:
        made += supp2(cr, a.out, a, order, labels, colours, cmap)
    if "3" in want:
        made += supp3(cr, a.out, a, order, labels, colours, cmap)
    if not made:
        die(f"--figures {a.figures!r} selected nothing; use 1, 2 and/or 3")

    spath = os.path.join(a.out, "supplemental_sources.tsv")
    with open(spath, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["panel", "source", "rows_used", "detail"])
        for row in SOURCES:
            w.writerow(row)
        w.writerow(["_run_", "jitter seed", a.seed, f"harmonic window +/-{HARMONIC_WINDOW} bp"])
        for m in WARNINGS:
            w.writerow(["_warning_", "", "", m])
    log(f"wrote {spath}")
    log(f"done: {len(made)} file(s); {len(WARNINGS)} warning(s)")


if __name__ == "__main__":
    main()
