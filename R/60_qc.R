#!/usr/bin/env Rscript
## Quality control: mappability and fragment (insert) sizes.
## Tables:  qc01_peak_mappability_models, qc01b_peak_mappability_by_sample, qc02_censat_mappability_stats,
##          qc03_insert_size_stats
## Figures: qc01_peak_mappability_vs_score, qc02_peak_mappability_by_censat, qc03_censat_base_mappability,
##          qc04_insert_<sample> (one per IP sample), qc05_insert_median_shift
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
cfg <- args_config()
qcc <- cfg$qc %||% list()
if (isFALSE(qcc$enabled)) { msg("qc disabled"); quit(save = "no") }
k <- qcc$mappability_k %||% 150
Q <- file.path(cfg$work, "qc")
S <- samples_dt(cfg)
clear_figures(cfg, "^qc")
cls_name <- function(x) ifelse(x == "outside_censat", "Outside CenSat regions",
                               ifelse(x == "genome", "Whole genome", pretty_censat(x)))
wrap <- function(x, w = 26) stringr_wrap(x, w)

## =============================== 1. peak mappability vs peak score / CenSat class ==============
pm <- if (file.exists(file.path(Q, "peak_mappability.tsv.gz"))) fread(cmd = paste("gzip -dc", shQuote(file.path(Q, "peak_mappability.tsv.gz")))) else NULL
if (!is.null(pm) && nrow(pm)) {
  ## peak scores from the original MACS3 files (column chosen by qc.peak_score)
  col <- switch(qcc$peak_score %||% "signal", signal = 7L, pvalue = 8L, qvalue = 9L, score = 5L,
                stop("qc.peak_score must be signal, pvalue, qvalue or score"))
  score_name <- switch(qcc$peak_score %||% "signal", signal = "fold enrichment", pvalue = "-log10 p-value",
                       qvalue = "-log10 q-value", score = "MACS3 score")
  sc <- rbindlist(lapply(seq_len(nrow(S)), function(i) {
    f <- S$peaks_file[i]
    if (!nzchar(f) || !file.exists(f)) return(NULL)
    d <- fread(f, header = FALSE, select = c(4L, col), col.names = c("name", "peak_score"))
    d[, peak_id := paste0(S$sample_id[i], "|", name)][, .(peak_id, peak_score)]
  }))
  pm <- merge(pm, sc, by = "peak_id", all.x = TRUE)
  ## CenSat class of each peak: reuse the labels from the main repeat analysis (step 06)
  pl <- rbindlist(lapply(list.files(file.path(cfg$work, "repeats"), "\\.censat\\.peak_labels\\.tsv$", full.names = TRUE),
                         function(f) fread(f, header = FALSE, col.names = c("peak_id", "chrom", "label"))))
  if (nrow(pl)) pm <- merge(pm, unique(pl[, .(peak_id, label)], by = "peak_id"), by = "peak_id", all.x = TRUE) else pm[, label := NA_character_]
  pm[is.na(label), label := "outside_censat"]
  pm[, category := cls_name(label)]
  pm <- pm[is.finite(peak_score) & peak_score > 0]
  pm[, log2_score := log2(peak_score)]
  lev <- c("Outside CenSat regions", setdiff(pm[, .N, by = category][order(-N), category], "Outside CenSat regions"))
  pm[, category := factor(category, levels = lev)]
  pm[, sample_id := factor(sample_id)]

  fit <- function(form, name) {
    m <- lm(form, data = pm)
    co <- summary(m)$coefficients
    d <- data.table(model = name, term = rownames(co), estimate = co[, 1], std_error = co[, 2], p_value = co[, 4])
    d <- d[!grepl("^\\(Intercept\\)|^sample_id", term)]
    d[, term_label := fifelse(term == "log2_score", sprintf("Peak score (%s), per doubling", score_name),
                       sub("^category", "CenSat class: ", paste0(term, " vs outside CenSat")))]
    d[, `:=`(ci_low = estimate - 1.96 * std_error, ci_high = estimate + 1.96 * std_error,
             n_peaks = nobs(m), r_squared = summary(m)$r.squared)]
    d
  }
  multi <- uniqueN(pm$sample_id) > 1
  adj <- if (multi) " + sample_id" else ""
  has_score <- uniqueN(pm$log2_score) > 1
  has_cat <- uniqueN(pm$category) > 1
  if (!has_score) msg("NOTE: all peaks have the same score; models and figure use CenSat class only")
  mods <- rbind(if (has_score) fit(as.formula(paste("mean_mappability ~ log2_score", adj)), "mappability ~ peak score"),
                if (has_cat) fit(as.formula(paste("mean_mappability ~ category", adj)), "mappability ~ CenSat class"),
                if (has_score && has_cat) fit(as.formula(paste("mean_mappability ~ log2_score + category", adj)), "mappability ~ peak score + CenSat class"))
  mods[, q_value := bh(p_value), by = model]
  mods[, note := if (multi) "sample included as a fixed effect (absorbs between-sample differences in score scale)" else ""]
  setcolorder(mods, c("model", "term", "term_label", "estimate", "std_error", "ci_low", "ci_high", "p_value", "q_value", "n_peaks", "r_squared", "note"))
  wt(mods, cfg, "qc01_peak_mappability_models")

  bys <- pm[, {
    ct <- if (.N >= 5 && uniqueN(log2_score) > 1) suppressWarnings(cor.test(mean_mappability, log2_score, method = "spearman", exact = FALSE)) else NULL
    .(n_peaks = .N, median_mappability = median(mean_mappability), pct_peaks_mappability_below_0.5 = 100 * mean(mean_mappability < 0.5),
      pct_peaks_in_censat = 100 * mean(category != "Outside CenSat regions"),
      spearman_rho_score = if (is.null(ct)) NA_real_ else unname(ct$estimate), p_value = if (is.null(ct)) NA_real_ else ct$p.value)
  }, by = sample_id]
  bys[, q_value := bh(p_value)]
  wt(bys, cfg, "qc01b_peak_mappability_by_sample")

  topc <- head(levels(pm$category), 12)
  pp <- pm[category %in% topc]
  if (has_score) {
  pp[, panel := factor(wrap(as.character(category)), levels = wrap(topc))]
  p1 <- ggplot(pp, aes(log2_score, mean_mappability)) +
    geom_bin2d(bins = 50) + scale_fill_gradient(low = "#e6f1fb", high = "#0c447c", trans = "log10", name = "Peaks") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#e34948", se = FALSE, linewidth = 0.7) +
    facet_wrap(~panel, ncol = 4) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = sprintf("Peak score: log2 %s", score_name), y = sprintf("Peak mean per-base mappability (k = %d)", k),
         title = "Peak mappability versus peak score, by CenSat class") + theme_x(9) + theme(legend.position = "right")
  save_plot(p1, cfg, "qc01_peak_mappability_vs_score", 13, 1.5 + 2.6 * ceiling(length(topc) / 4))
  b <- mods[model == "mappability ~ peak score" & term == "log2_score"]
  fig_note(cfg, "qc01_peak_mappability_vs_score",
           paste0(sprintf("%d peaks from %d samples; %d CenSat classes with the most peaks shown.", nrow(pm), uniqueN(pm$sample_id), length(topc)),
                  if (nrow(b) == 1) sprintf(" Overall (model with sample as a fixed effect): mappability changes by %.3f per doubling of peak score (95%% CI %.3f to %.3f).",
                                            b$estimate, b$ci_low, b$ci_high) else ""))
  } else fig_note(cfg, "qc01__summary", "Not drawn: every peak has the same score, so mappability cannot be related to it.")

  cnt <- pm[, .N, by = category]
  p2 <- ggplot(pm, aes(category, mean_mappability)) +
    geom_boxplot(fill = "#b5d4f4", outlier.size = 0.3, outlier.alpha = 0.3) +
    geom_text(data = cnt, aes(y = 1.04, label = paste0("n=", N)), size = 2.8) +
    scale_x_discrete(labels = function(x) wrap(x, 18)) + coord_cartesian(ylim = c(0, 1.06)) +
    labs(x = NULL, y = sprintf("Peak mean per-base mappability (k = %d)", k), title = "Peak mappability by CenSat class (all samples)") +
    theme_x(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_plot(p2, cfg, "qc02_peak_mappability_by_censat", max(9, 2 + 0.7 * nrow(cnt)), 6)
  fig_note(cfg, "qc02_peak_mappability_by_censat", sprintf("%d peaks from %d samples.", nrow(pm), uniqueN(pm$sample_id)))
}

## =============================== 2. per-base mappability statistics per CenSat class ===========
h <- rd(file.path(Q, "censat_mappability_hist.tsv"))
if (!is.null(h) && nrow(h)) {
  h[, m := mappable_kmers / k]
  st <- h[, {
    o <- order(m); cw <- cumsum(n_bases[o]) / sum(n_bases)
    .(total_bp = sum(n_bases), mean_mappability = sum(m * n_bases) / sum(n_bases),
      median_mappability = m[o][which(cw >= 0.5)[1]],
      pct_bases_fully_unique = 100 * sum(n_bases[mappable_kmers == k]) / sum(n_bases),
      pct_bases_mappability_ge_0.5 = 100 * sum(n_bases[m >= 0.5]) / sum(n_bases),
      pct_bases_unmappable = 100 * sum(n_bases[mappable_kmers == 0]) / sum(n_bases))
  }, by = class]
  st[, class_name := cls_name(class)]
  if (exists("pm") && !is.null(pm) && nrow(pm)) {
    ps <- pm[, .(n_peaks_all_samples = .N, median_peak_mappability = median(mean_mappability)), by = label]
    st <- merge(st, ps, by.x = "class", by.y = "label", all.x = TRUE)
  }
  setcolorder(st, c("class", "class_name"))
  setorder(st, mean_mappability)
  wt(st, cfg, "qc02_censat_mappability_stats")

  h[, bin := cut(m, breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 0.999999, Inf),
                 labels = c("0 (no single-copy k-mers)", "0-0.25", "0.25-0.5", "0.5-0.75", "0.75-<1", "1 (fully unique)"))]
  hb <- h[, .(bp = sum(n_bases)), by = .(class, bin)][, pct := 100 * bp / sum(bp), by = class]
  hb[, class_name := factor(cls_name(class), levels = cls_name(st$class))]
  p3 <- ggplot(hb, aes(class_name, pct, fill = bin)) + geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("#a32d2d", "#eb6834", "#eda100", "#b5d4f4", "#5a9be0", "#0c447c"), name = "Per-base mappability", drop = FALSE) +
    scale_x_discrete(labels = function(x) wrap(x, 40)) + guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = "% of bases", title = sprintf("Per-base mappability of CenSat classes (k = %d, single-copy k-mers)", k)) + theme_x(10)
  save_plot(p3, cfg, "qc03_censat_base_mappability", 11, 2.5 + 0.35 * uniqueN(hb$class))
  fig_note(cfg, "qc03_censat_base_mappability", sprintf("All bases of each class genome-wide; classes ordered by mean mappability. Genome-wide, %.1f%% of bases are fully unique.",
                                                       st[class == "genome", pct_bases_fully_unique]))
}

## =============================== 3. insert sizes per CenSat class, IP vs control ===============
libs <- rd(file.path(Q, "insert", "libraries.tsv"), header = FALSE, col.names = c("library_id", "bam"))
if (!is.null(libs) && nrow(libs)) {
  ins <- rbindlist(lapply(libs$library_id, function(id) {
    d <- rd(file.path(Q, "insert", paste0(id, ".tsv")))
    if (is.null(d) || !nrow(d)) return(NULL)
    d[, library_id := id]
  }))
  status <- vapply(libs$library_id, function(id) { f <- file.path(Q, "insert", paste0(id, ".status")); if (file.exists(f)) readLines(f)[1] else "missing" }, "")
  if (nrow(ins)) {
    ins[, class_name := cls_name(class)]
    roles <- rbind(S[bam_file != "", .(library_id = libs$library_id[match(bam_file, libs$bam)], role = "IP", sample_id)],
                   S[control_bam != "", .(library_id = libs$library_id[match(control_bam, libs$bam)], role = "control", sample_id)])
    used <- roles[, .(role = paste(unique(role), collapse = "/"), used_by = paste(unique(sample_id), collapse = ",")), by = library_id]
    stt <- ins[, {
      o <- order(insert_size); cw <- cumsum(n[o]) / sum(n); x <- insert_size[o]
      .(n_fragments = sum(n), mean_insert = sum(insert_size * n) / sum(n), q25 = x[which(cw >= .25)[1]],
        median_insert = x[which(cw >= .5)[1]], q75 = x[which(cw >= .75)[1]], mode_insert = insert_size[which.max(n)])
    }, by = .(library_id, class, class_name)]
    stt <- merge(used, stt, by = "library_id")
    setorder(stt, library_id, -n_fragments)
    wt(stt, cfg, "qc03_insert_size_stats")

    ## one figure per IP sample: IP and its control, per class (top 11 classes by IP fragments + outside)
    pairs <- S[bam_file != "", .(sample_id, ip = libs$library_id[match(bam_file, libs$bam)], ctl = libs$library_id[match(control_bam, libs$bam)])]
    no_ctl <- character()
    for (i in seq_len(nrow(pairs))) {
      s <- pairs$sample_id[i]; ip <- pairs$ip[i]; ctl <- pairs$ctl[i]
      d <- ins[library_id %in% c(ip, ctl)]
      if (!nrow(d[library_id == ip])) next
      if (is.na(ctl) || !nrow(d[library_id == ctl])) no_ctl <- c(no_ctl, s)
      d[, library := fifelse(library_id == ip, "IP", "Control")]
      topc <- c(head(setdiff(d[library == "IP", .(t = sum(n)), by = class][order(-t), class], "outside_censat"), 11), "outside_censat")
      d <- d[class %in% topc]
      d[, bin := 10 * (insert_size %/% 10)]
      d <- d[, .(n = sum(n)), by = .(library, class, class_name, bin)][, frac := n / sum(n), by = .(library, class)]
      totals <- d[, .(t = sum(n)), by = .(library, class)]
      d[, panel := factor(wrap(paste0(class_name, "\n(IP ", format(totals[library == "IP" & class == .BY$class, t], big.mark = ","),
                                    " / control ", format(totals[library == "Control" & class == .BY$class, t][1], big.mark = ","), " fragments)"), 40)),
        by = class]
      d[, panel := factor(panel, levels = unique(panel[order(match(class, topc))]))]
      p4 <- ggplot(d, aes(bin, frac, colour = library)) + geom_line(linewidth = 0.7) +
        scale_colour_manual(values = c(IP = "#e34948", Control = "#2a78d6"), name = NULL) +
        facet_wrap(~panel, ncol = 4, scales = "free_y") +
        labs(x = "Insert (fragment) size, bp (10-bp bins)", y = "Fraction of fragments in the class",
             title = sprintf("Insert-size distributions by CenSat class: %s", s)) + theme_x(8)
      fn <- paste0("qc04_insert_", gsub("[^A-Za-z0-9]+", "_", s))
      save_plot(p4, cfg, fn, 13, 1.5 + 2.4 * ceiling(length(topc) / 4))
      fig_note(cfg, fn, sprintf("IP library %s versus control library %s; %d classes shown (most IP fragments first, plus all fragments outside CenSat, subsampled).",
                                ip, ifelse(is.na(ctl), "none", ctl), length(topc)))
    }
    se <- names(status)[status == "single_end"]
    fig_note(cfg, "qc04__summary", paste(c(if (length(se)) sprintf("Single-end libraries (no insert sizes): %s.", paste(se, collapse = ", ")),
                                           if (length(no_ctl)) sprintf("No paired-end control for: %s.", paste(no_ctl, collapse = ", "))), collapse = " "))

    ## summary: median IP insert minus median control insert, per class and sample
    med <- stt[, .(library_id, class, class_name, median_insert, n_fragments)]
    sh <- merge(pairs, med, by.x = "ip", by.y = "library_id", allow.cartesian = TRUE)
    sh <- merge(sh, med[, .(ctl = library_id, class, median_ctl = median_insert, n_ctl = n_fragments)], by = c("ctl", "class"))
    sh <- sh[n_fragments >= 50 & n_ctl >= 50]
    if (nrow(sh)) {
      sh[, shift := median_insert - median_ctl]
      cl_order <- sh[, .(t = sum(n_fragments)), by = class_name][order(t), class_name]
      p5 <- ggplot(sh, aes(short_id(sample_id), factor(class_name, levels = cl_order), fill = shift)) + geom_tile(colour = "white") +
        geom_text(aes(label = round(shift)), size = 2.4) +
        scale_fill_gradient2(low = "#2a78d6", mid = "#f7f6f2", high = "#e34948", name = "IP - control\nmedian (bp)") +
        labs(x = NULL, y = NULL, title = "Median insert size, IP minus control, by CenSat class") +
        theme_x(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
      save_plot(p5, cfg, "qc05_insert_median_shift", max(9, 4 + 0.35 * uniqueN(sh$sample_id)), 2.5 + 0.35 * length(cl_order))
      fig_note(cfg, "qc05_insert_median_shift", "Cells with fewer than 50 fragments in either library are left blank.")
    }
  }
}
