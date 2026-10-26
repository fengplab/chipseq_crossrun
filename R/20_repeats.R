#!/usr/bin/env Rscript
## Q2: Do peaks fall on known repeats / centromeric satellites more than expected?
## Tables:  table03_repeat_enrichment (long), table03_log2oe_<featureset> (wide),
##          table03b_censat_by_chromosome, table03c_repeat_enrichment_vs_universe,
##          table03d_reads_in_repeats (if step 07 ran)
## Figures: fig04_log2oe_<featureset>_heatmap, fig05_censat_composition, fig06_censat_by_chromosome,
##          fig07_reads_in_repeats_heatmap
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
suppressPackageStartupMessages(library(pheatmap))
cfg <- args_config()
rc <- cfg$repeats
pc <- rc$pseudocount %||% 1
S <- samples_dt(cfg)
ps <- rd(file.path(cfg$work, "peaks", "peak_stats.tsv"))

cnt <- rbindlist(lapply(list.files(file.path(cfg$work, "repeats"), "\\.counts\\.tsv$", full.names = TRUE), rd))
if (!nrow(cnt)) stop("no repeat counts; run step 06")
n_iter <- cnt[iter > 0, uniqueN(iter)]

## every label x sample x iteration must exist (absent = 0)
grid <- CJ(sample_id = unique(cnt$sample_id), iter = 0:n_iter, sorted = FALSE)
labs <- unique(cnt[, .(feature_set, label)])
full <- grid[, as.list(labs), by = .(sample_id, iter)]
full <- merge(full, cnt, by = c("sample_id", "iter", "feature_set", "label"), all.x = TRUE)
full[is.na(n), n := 0]

obs <- full[iter == 0, .(sample_id, feature_set, label, observed = n)]
nul <- full[iter > 0, .(expected = mean(n), expected_sd = sd(n)), by = .(sample_id, feature_set, label)]
ge <- merge(full[iter > 0], obs, by = c("sample_id", "feature_set", "label"))[
  , .(p_enrich = (sum(n >= observed) + 1) / (.N + 1), p_deplete = (sum(n <= observed) + 1) / (.N + 1)),
  by = .(sample_id, feature_set, label)]
t3 <- Reduce(function(a, b) merge(a, b, by = c("sample_id", "feature_set", "label")), list(obs, nul, ge))
t3 <- merge(t3, ps[, .(sample_id, total_peaks = n_peaks)], by = "sample_id", all.x = TRUE)
t3[, `:=`(pct_of_peaks = 100 * observed / pmax(total_peaks, 1),
          log2_oe = log2((observed + pc) / (expected + pc)),
          z = fifelse(expected_sd > 0, (observed - expected) / expected_sd, NA_real_))]
t3[, `:=`(q_enrich = bh(p_enrich), q_deplete = bh(p_deplete)), by = feature_set]
settings <- rd(file.path(cfg$work, "repeats", "null_settings.tsv"), header = FALSE)
if (!is.null(settings)) for (i in seq_len(nrow(settings))) t3[, (paste0("null_", settings$V1[i])) := settings$V2[i]]
t3[, label_name := pretty_label(feature_set, label)]
setorder(t3, feature_set, label, sample_id)
wt(t3, cfg, "table03_repeat_enrichment")

min_obs <- rc$min_observed %||% 5
for (fs in unique(t3$feature_set)) {
  d <- t3[feature_set == fs]
  ## rows: for every sample, the classes significantly enriched or depleted (q < 0.05), or the 5 with the smallest
  ## p-values when fewer are significant; union over samples (classes need >= min_observed peaks in some sample)
  cand <- d[label != "ANY", .(mx = max(observed)), by = label][mx >= min_obs, label]
  dd <- d[label %in% cand][, `:=`(p_any = pmin(p_enrich, p_deplete), is_sig = q_enrich < 0.05 | q_deplete < 0.05)]
  keep <- select_heatmap_rows(dd, "label", "sample_id", "p_any", "is_sig", 5)
  w <- dcast(d[label %in% c(keep, "ANY")], label ~ sample_id, value.var = "log2_oe")
  wt(w, cfg, paste0("table03_log2oe_", fs))
  m <- as.matrix(w[label != "ANY", -1]); rownames(m) <- w[label != "ANY", label]
  if (nrow(m) < 1 || ncol(m) < 1) next
  raw_rows <- rownames(m)
  lim <- max(1, min(6, max(abs(m), na.rm = TRUE)))
  d[, q_any := pmin(q_enrich, q_deplete)]
  qs <- dcast(d[label %in% rownames(m)], label ~ sample_id, value.var = "q_any")
  star <- as.matrix(qs[match(rownames(m), qs$label), -1])[, colnames(m), drop = FALSE]
  star <- ifelse(!is.na(star) & star < 0.05, "*", "")
  rownames(m) <- pretty_label(fs, raw_rows); rownames(star) <- rownames(m)
  fs_name <- if (fs == "censat") "Centromeric/satellite (CenSat) classes" else "RepeatMasker classes"
  ## Size the page from the text so nothing is clipped. pheatmap centres the title on the heatmap body and reserves
  ## one line for it, so the title stays on one line and the body is made at least as wide as the title; the rotated
  ## column labels (without the _REP suffix) and the row labels get the room they need.
  ann_c <- sample_annotation(cfg, S, colnames(m))
  ann_w <- if (is.null(ann_c)) 0 else max(text_in(names(ann_c), 8), text_in(unlist(lapply(ann_c, function(v) unique(as.character(v)))), 8)) + 0.45
  ttl <- sprintf("%s: peak enrichment, log2(observed / expected); * q < 0.05", if (fs == "censat") "CenSat classes" else "RepeatMasker classes")
  side_w <- text_in(rownames(m), 8) + ann_w + 1.6 + (if (nrow(m) > 2) 0.7 else 0)   # row labels + legends + colour key + margins (+ row tree)
  fig_w <- max(9, side_w + max(0.3 * ncol(m), text_in(ttl, 8 * 1.3) + 0.2))
  ## fixed cell sizes: pheatmap otherwise stretches the body over the page and the rotated labels run off the bottom
  body_w <- fig_w - side_w; body_h <- max(1.2, 0.22 * nrow(m))
  fig_h <- body_h + text_in(short_id(colnames(m)), 8) + 0.17 * length(ann_c) + 1.5    # + title, dendrogram, margins
  save_pheatmap(function() pheatmap(m, cluster_rows = nrow(m) > 2, cluster_cols = ncol(m) > 2,
                                    breaks = seq(-lim, lim, length.out = 51),
                                    color = colorRampPalette(c("#2a78d6", "#f0efec", "#e34948"))(50),
                                    display_numbers = star, annotation_col = ann_c, labels_col = short_id(colnames(m)),
                                    cellwidth = body_w * 72 / ncol(m), cellheight = body_h * 72 / nrow(m),
                                    fontsize = 8, main = ttl),
                cfg, paste0("fig04_log2oe_", fs, "_heatmap"), width = fig_w, height = fig_h)
  fig_note(cfg, paste0("fig04_log2oe_", fs, "_heatmap"),
           sprintf("%d classes: for each sample, every class enriched or depleted with q < 0.05, or the 5 with the smallest p-values when fewer are significant (classes with >= %d peaks in at least one sample).", nrow(m), min_obs))
}

## ---------------- CenSat composition and chromosome drill-down --------------------------------------
cs <- t3[feature_set == "censat" & label != "ANY"]
if (nrow(cs)) {
  topl <- cs[, .(tot = sum(observed)), by = label][order(-tot)][seq_len(min(8, .N)), label]
  cs[, lab := fifelse(label %in% topl, pretty_censat(label), "Other CenSat classes")]
  cc <- cs[, .(pct = sum(pct_of_peaks)), by = .(sample_id, lab)]
  lev <- c(pretty_censat(topl), "Other CenSat classes")
  cc[, lab := factor(lab, levels = lev)]
  pal <- setNames(c(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")[seq_along(topl)], "grey70"), lev)
  ## Size from the text: 45-degree sample labels hang left of their bar (the first one past the panel edge unless the
  ## left margin makes room) and need height below the axis; the legend needs its full width, else it is cut at the sides.
  xl <- short_id(unique(cc$sample_id)); ax_pt <- 8.8                                # theme_minimal(10) axis text
  lab_run <- text_in(xl, ax_pt) * cos(pi / 4)                                       # horizontal = vertical extent of a label
  key_w <- text_in(lev, 8.8) + 0.45                                                 # one legend column: key + text
  w5 <- max(8, 2 + 0.35 * length(xl))
  ncol5 <- if (2 * key_w + 0.4 <= w5) 2 else 1
  w5 <- max(w5, ncol5 * key_w + 0.6)
  h5 <- 4.6 + lab_run + 0.22 * ceiling(length(lev) / ncol5)
  first_run <- text_in(xl[order(xl)][1], ax_pt) * cos(pi / 4) * 72                  # pt; x is sorted alphabetically
  p5 <- ggplot(cc, aes(short_id(sample_id), pct, fill = lab)) + geom_col(width = 0.7) +
    scale_fill_manual(values = pal) + guides(fill = guide_legend(ncol = ncol5)) +
    labs(x = NULL, y = "% of peaks", fill = NULL, title = "Peaks in centromeric / satellite (CenSat) classes") +
    theme_x() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      plot.margin = margin(6, 12, 6, max(6, first_run - 0.5 * (w5 * 72 - 60) / length(xl))))
  save_plot(p5, cfg, "fig05_censat_composition", w5, h5)
  fig_note(cfg, "fig05_censat_composition", sprintf("The %d CenSat classes with the most peaks overall are shown individually; the rest are pooled as 'Other'.", length(topl)))

  pl <- rbindlist(lapply(list.files(file.path(cfg$work, "repeats"), "\\.censat\\.peak_labels\\.tsv$", full.names = TRUE),
                         function(f) { d <- fread(f, header = FALSE, col.names = c("peak_id", "chrom", "label"))
                                       if (nrow(d)) d[, sample_id := sub("\\|.*$", "", peak_id)]; d }), fill = TRUE)
  if (nrow(pl)) {
    ## totals across all samples (no per-sample breakdown)
    bych <- pl[, .(n_peaks = .N, n_samples = uniqueN(sample_id)), by = .(chrom, label)]
    bych[, censat_class := pretty_censat(label)]
    bych[, pct_of_class_peaks := 100 * n_peaks / sum(n_peaks), by = label]
    chr_order <- unique(bych$chrom)[order(suppressWarnings(as.numeric(sub("^chr", "", unique(bych$chrom)))), unique(bych$chrom))]
    setorder(bych, label, -n_peaks)
    wt(bych[, .(chrom, label, censat_class, n_peaks, n_samples, pct_of_class_peaks)], cfg, "table03b_censat_by_chromosome")
    ## every chromosome of the genome, in genome order, including those without CenSat peaks (shown as 0)
    gs <- rd(file.path(cfg$work, "annot", "genome.sizes"), header = FALSE)
    all_chr <- if (!is.null(gs)) gs$V1 else unique(bych$chrom)
    chr_order <- all_chr[order(suppressWarnings(as.numeric(sub("^chr", "", all_chr))), all_chr)]
    cls_order <- bych[, .(t = sum(n_peaks)), by = censat_class][order(t), censat_class]
    bych <- merge(CJ(chrom = chr_order, censat_class = cls_order), bych[, .(chrom, censat_class, n_peaks)],
                  by = c("chrom", "censat_class"), all.x = TRUE)
    bych[is.na(n_peaks), n_peaks := 0]
    bych[, chrom := factor(chrom, levels = chr_order)]
    bych[, censat_class := factor(censat_class, levels = cls_order)]
    p6 <- ggplot(bych, aes(chrom, censat_class, fill = ifelse(n_peaks > 0, n_peaks, NA))) + geom_tile(colour = "white") +
      geom_text(aes(label = n_peaks), size = 2.4, colour = "grey15") +
      scale_fill_gradient(low = "#e6f1fb", high = "#2a78d6", trans = "log10", na.value = "grey95", name = "Peaks (all samples)") +
      labs(x = NULL, y = NULL, title = "Peaks in CenSat classes by chromosome, summed over all samples") + theme_x(9) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "right")
    save_plot(p6, cfg, "fig06_censat_by_chromosome", max(10, 3 + 0.45 * length(chr_order)), 2 + 0.35 * length(cls_order))
    fig_note(cfg, "fig06_censat_by_chromosome", sprintf("%d peaks in CenSat classes across %d samples; all %d chromosomes shown (grey = no peaks).", sum(bych$n_peaks), uniqueN(pl$sample_id), length(chr_order)))
  }
}

## ---------------- sample vs pooled universe (LOLA-style) ---------------------------------------------
## Is this sample more enriched than the typical ChIP in this collection? 2x2 per label:
## (this sample's peaks in / not in label) vs (all merged regions across samples in / not in label).
rl <- rd(file.path(cfg$work, "regions", "region_labels.tsv"))
mr <- rd(file.path(cfg$work, "regions", "merged_regions.tsv"))
if (!is.null(rl) && !is.null(mr)) {
  U <- nrow(mr)
  ul <- rl[, .(u_in = uniqueN(region_id)), by = .(feature_set, label)]
  uv <- merge(t3[label != "ANY", .(sample_id, feature_set, label, observed, total_peaks)], ul, by = c("feature_set", "label"))
  uv[, c("odds_ratio", "p_value") := {
    r <- mapply(function(a, n, b) {
      ft <- fisher.test(matrix(c(a, n - a, b, U - b), 2))
      c(unname(ft$estimate), ft$p.value)
    }, observed, total_peaks, u_in)
    list(r[1, ], r[2, ])
  }]
  uv[, `:=`(pct_sample = 100 * observed / total_peaks, pct_universe = 100 * u_in / U, q_value = bh(p_value))]
  wt(uv[order(q_value)], cfg, "table03c_repeat_enrichment_vs_universe")
}

## ---------------- read-level signal (optional step 07) ----------------------------------------------
rr <- rd(file.path(cfg$work, "reads", "reads_in_labels.tsv"))
if (!is.null(rr) && nrow(rr)) {
  rr[, frac := reads / total_reads]
  ip <- S[bam_file != "", .(sample_id, bam = bam_file, control_bam)]
  ipm <- merge(ip, rr, by = "bam")
  ctl <- rr[, .(control_bam = bam, feature_set, label, control_frac = frac)]
  ipm <- merge(ipm, ctl, by = c("control_bam", "feature_set", "label"), all.x = TRUE)
  ipm[, log2_ip_over_control := log2((frac + 1e-7) / (control_frac + 1e-7))]
  wt(ipm[, .(sample_id, feature_set, label, ip_reads = reads, ip_total = total_reads, ip_frac = frac,
             control_frac, log2_ip_over_control)], cfg, "table03d_reads_in_repeats")
  for (fs in unique(ipm$feature_set)) {
    dd <- ipm[feature_set == fs]
    dd[!is.finite(log2_ip_over_control), log2_ip_over_control := NA_real_]
    if (fs == "censat") {
      ## every sample and every CenSat class, also those without reads or without a control (blank cells)
      cls <- unique(c(dd$label, tryCatch(fread(file.path(cfg$work, "annot", "censat.class.bed"), header = FALSE, select = 4)[[1]], error = function(e) character())))
      dd <- merge(CJ(label = setdiff(cls, "ANY"), sample_id = S$sample_id), dd[, .(label, sample_id, log2_ip_over_control)],
                  by = c("label", "sample_id"), all.x = TRUE)
    }
    w <- dcast(dd, label ~ sample_id, value.var = "log2_ip_over_control")
    m <- as.matrix(w[, -1]); rownames(m) <- pretty_label(fs, w$label)
    if (fs != "censat") m <- m[rowSums(is.na(m)) < ncol(m), , drop = FALSE]
    if (nrow(m) < 2 || ncol(m) < 1) next
    has_na <- anyNA(m)   # pheatmap cannot cluster with missing cells; keep a fixed order then
    fs_name <- if (fs == "censat") "CenSat classes" else "RepeatMasker classes"
    lim <- max(1, min(5, max(abs(m), na.rm = TRUE)))
    if (has_na) m <- m[order(rownames(m)), , drop = FALSE]
    save_pheatmap(function() pheatmap(m, cluster_rows = !has_na, cluster_cols = !has_na && ncol(m) > 2, breaks = seq(-lim, lim, length.out = 51), na_col = "grey90",
                                      color = colorRampPalette(c("#2a78d6", "#f0efec", "#e34948"))(50), fontsize = 8,
                                      annotation_col = sample_annotation(cfg, S, colnames(m)),
                                      main = sprintf("%s: log2(IP read fraction / control read fraction)", fs_name)),
                  cfg, paste0("fig07_reads_in_", fs, "_heatmap"), width = max(9, 4 + 0.3 * ncol(m)), height = max(5, 3 + 0.22 * nrow(m)))
    fig_note(cfg, paste0("fig07_reads_in_", fs, "_heatmap"), sprintf("%d classes x %d samples%s.", nrow(m), ncol(m),
             if (has_na) "; grey = no reads in the class or no control library (rows and columns are not clustered when cells are missing)" else ""))
  }
}
