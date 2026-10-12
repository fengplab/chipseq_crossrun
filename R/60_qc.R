#!/usr/bin/env Rscript
## Quality control: mappability and fragment (insert) sizes.
## Tables:  qc01_peak_mappability_models, qc01b_peak_mappability_by_sample, qc01c_mappability_score_by_censat,
##          qc02_censat_mappability_stats, qc03_insert_size_stats
## Figures: qc01_peak_mappability_vs_score, qc02_peak_mappability_by_censat (violins),
##          qc03_insert_<sample> (one per IP sample), qc04_insert_median_shift
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
  if (!nrow(sc)) {
    miss <- S[nzchar(peaks_file) & !file.exists(peaks_file), peaks_file]
    stop("none of the MACS3 peak files listed in work/samples.tsv can be read",
         if (length(miss)) paste0(" (e.g. ", miss[1], "). If the pipeline results moved, update results_dir in the manifest and rerun step 0.") else "",
         call. = FALSE)
  }
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
  ## per-class linear trend of mappability on log2 score; BH across classes; line drawn only if q < 0.05
  tr <- pp[, {
    if (.N >= 10 && uniqueN(log2_score) > 2 && uniqueN(mean_mappability) > 2) {
      co <- summary(lm(mean_mappability ~ log2_score))$coefficients
      .(n_peaks = .N, slope_per_doubling = co[2, 1], std_error = co[2, 2], p_value = co[2, 4])
    } else .(n_peaks = .N, slope_per_doubling = NA_real_, std_error = NA_real_, p_value = NA_real_)
  }, by = category]
  tr[, q_value := bh(p_value)]
  tr[, line_drawn := !is.na(q_value) & q_value < 0.05]
  wt(tr[, .(censat_class = as.character(category), n_peaks, slope_per_doubling, std_error, p_value, q_value, line_drawn)][order(q_value)],
     cfg, "qc01c_mappability_score_by_censat")
  pp[, panel := factor(wrap(as.character(category)), levels = wrap(topc))]
  tr[, panel := factor(wrap(as.character(category)), levels = wrap(topc))]
  tr[, lab := fifelse(is.na(p_value), fifelse(n_peaks < 10, "too few peaks to test", "no variation to test"),
                      sprintf("slope %+.2g per doubling, q %s", slope_per_doubling, ifelse(q_value < 2.2e-16, fmt_q(q_value), paste("=", fmt_q(q_value)))))]
  sig_cats <- tr[line_drawn == TRUE, category]
  p1 <- ggplot(pp, aes(log2_score, mean_mappability)) +
    geom_bin2d(bins = 50) + scale_fill_gradient(low = "#e6f1fb", high = "#0c447c", trans = "log10", name = "Peaks") +
    geom_smooth(data = pp[category %in% sig_cats], method = "lm", formula = y ~ x, colour = "#e34948", se = FALSE, linewidth = 0.8) +
    geom_text(data = tr, aes(label = lab), x = -Inf, y = Inf, hjust = -0.03, vjust = 1.4, size = 2.6, inherit.aes = FALSE) +
    facet_wrap(~panel, ncol = 4) + coord_cartesian(ylim = c(0, 1.08)) +
    labs(x = sprintf("Peak score: log2 %s", score_name), y = sprintf("Peak mean per-base mappability (k = %d)", k),
         title = "Peak mappability versus peak score, by CenSat class",
         subtitle = "Red line: linear trend, drawn only where significant (q < 0.05, BH across classes)") +
    theme_x(9) + theme(legend.position = "right")
  save_plot(p1, cfg, "qc01_peak_mappability_vs_score", 13, 1.5 + 2.6 * ceiling(length(topc) / 4))
  fig_note(cfg, "qc01_peak_mappability_vs_score",
           sprintf("%d peaks from %d samples; %d CenSat classes with the most peaks shown. A significant trend (q < 0.05) was found in %d of %d tested classes.",
                   nrow(pm), uniqueN(pm$sample_id), length(topc), length(sig_cats), sum(!is.na(tr$p_value))))
  } else fig_note(cfg, "qc01__summary", "Not drawn: every peak has the same score, so mappability cannot be related to it.")

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

}

## ---- QC2: violins of peak mappability per CenSat class, with the class's genome-wide median per-base
##      mappability (red line) and the peaks' own median (black dot)
if (exists("pm") && !is.null(pm) && nrow(pm)) {
  cnt <- pm[, .N, by = category]
  ref <- if (exists("st") && !is.null(h) && nrow(h)) st[, .(category = class_name, class_median = median_mappability)] else NULL
  if (!is.null(ref)) ref <- ref[category %in% levels(pm$category)][, category := factor(category, levels = levels(pm$category))]
  p2 <- ggplot(pm, aes(category, mean_mappability)) +
    geom_violin(data = pm[category %in% cnt[N >= 2, category]], fill = "#b5d4f4", colour = "#2a78d6", scale = "width", linewidth = 0.3) +
    stat_summary(fun = median, geom = "point", size = 1.8, colour = "black") +
    { if (!is.null(ref) && nrow(ref)) geom_errorbar(data = ref, aes(x = category, ymin = class_median, ymax = class_median),
                                                     inherit.aes = FALSE, colour = "#e34948", linewidth = 0.9, width = 0.85) } +
    geom_text(data = cnt, aes(y = 1.06, label = paste0("n=", N)), size = 2.8) +
    scale_x_discrete(labels = function(x) wrap(x, 18)) + coord_cartesian(ylim = c(0, 1.08)) +
    labs(x = NULL, y = sprintf("Peak mean per-base mappability (k = %d)", k), title = "Peak mappability by CenSat class (all samples)",
         subtitle = "Red line: median per-base mappability of the whole class genome-wide; black dot: median of the class's peaks") +
    theme_x(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_plot(p2, cfg, "qc02_peak_mappability_by_censat", max(10, 2 + 0.75 * nrow(cnt)), 6.5)
  fig_note(cfg, "qc02_peak_mappability_by_censat", sprintf("%d peaks from %d samples.%s", nrow(pm), uniqueN(pm$sample_id),
           if (is.null(ref)) " Genome-wide class medians unavailable (no mappability histogram)." else ""))
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

    ## one figure per IP sample: IP and its control, per class (top 11 classes by IP fragments + outside).
    ## These are extra figures (extra_plots.insert_sizes), written to <outdir>/extra/insert_sizes/ and not part of the report.
    clear_figures(cfg, "^qc03_insert_")
    fig_note(cfg, "qc03__summary", "")
    xdir <- file.path(cfg$outdir, "extra", "insert_sizes")
    do_ins <- isTRUE((cfg$extra_plots %||% list())$insert_sizes %||% TRUE)
    if (do_ins) dir.create(xdir, showWarnings = FALSE, recursive = TRUE)
    pairs <- S[bam_file != "", .(sample_id, ip = libs$library_id[match(bam_file, libs$bam)], ctl = libs$library_id[match(control_bam, libs$bam)])]
    no_ctl <- character()
    for (i in seq_len(nrow(pairs))) {
      if (!do_ins) break
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
      p4 <- p4 + labs(subtitle = sprintf("IP library %s versus control library %s; %d classes (most IP fragments first, plus fragments outside CenSat, subsampled)",
                                         ip, ifelse(is.na(ctl), "none", ctl), length(topc)))
      fn <- file.path(xdir, paste0(gsub("[^A-Za-z0-9]+", "_", s), ".png"))
      ggsave(fn, p4, width = 13, height = 1.5 + 2.4 * ceiling(length(topc) / 4), dpi = 150, bg = "white")
    }
    if (do_ins) msg(sprintf("insert-size figures for %d samples in %s", nrow(pairs), xdir))
    se <- names(status)[status == "single_end"]
    if (length(se)) msg("single-end libraries (no insert sizes): ", paste(se, collapse = ", "))
    if (length(no_ctl)) msg("no paired-end control for: ", paste(no_ctl, collapse = ", "))

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
      save_plot(p5, cfg, "qc04_insert_median_shift", max(9, 4 + 0.35 * uniqueN(sh$sample_id)), 2.5 + 0.35 * length(cl_order))
      fig_note(cfg, "qc04_insert_median_shift", "Cells with fewer than 50 fragments in either library are left blank.")
    }
  }
}

## =============================== 4. strand balance and read stacking in peaks (QC5) ============
## Real ChIP peaks draw fragments from both strands and from many distinct positions; alignment artifacts in
## repeats tend to be strand-biased and stacked on a few start positions.
shape_files <- list.files(file.path(Q, "peak_shape"), "\\.tsv$", full.names = TRUE)
if (length(shape_files)) {
  shp <- rbindlist(lapply(shape_files, rd))
  bed <- rbindlist(lapply(list.files(file.path(cfg$work, "peaks"), "\\.bed5$", full.names = TRUE),
                          function(f) fread(f, header = FALSE, select = 2:4, col.names = c("start", "end", "peak_id"))))
  shp <- merge(shp, bed, by = "peak_id")
  pl5 <- rbindlist(lapply(list.files(file.path(cfg$work, "repeats"), "\\.censat\\.peak_labels\\.tsv$", full.names = TRUE),
                          function(f) fread(f, header = FALSE, col.names = c("peak_id", "chrom", "label"))))
  shp <- merge(shp, if (nrow(pl5)) unique(pl5[, .(peak_id, label)], by = "peak_id") else data.table(peak_id = character(), label = character()),
               by = "peak_id", all.x = TRUE)
  shp[is.na(label), label := "outside_censat"]
  minr <- qcc$shape_min_reads %||% 10
  n_all <- nrow(shp)
  shp <- shp[n_reads >= minr]
  if (nrow(shp)) {
    shp[, width := pmax(1, end - start)]
    shp[, strand_imbalance := abs(n_forward / n_reads - 0.5) * 2]
    shp[, strand_p := pmin(1, 2 * pmin(pbinom(n_forward, n_reads, 0.5), pbinom(n_forward - 1, n_reads, 0.5, lower.tail = FALSE)))]
    shp[, M := 2 * width][, expected_distinct := M * (1 - (1 - 1 / M)^n_reads)]
    shp[, start_diversity := n_distinct_starts / expected_distinct]
    shp[, top_start_share := n_top_start / n_reads]
    shp[, category := cls_name(label)]
    lev5 <- c("Outside CenSat regions", setdiff(shp[, .N, by = category][order(-N), category], "Outside CenSat regions"))
    shp[, category := factor(category, levels = lev5)]
    mets <- c(strand_imbalance = "Strand imbalance: |forward share - 0.5| x 2 (0 = balanced, 1 = one strand)",
              start_diversity = "Start-site diversity: distinct starts / expected if reads were placed at random (1 = random, low = stacked)",
              top_start_share = "Share of the peak's reads at its most-used start position")
    ref <- shp[category == "Outside CenSat regions"]
    s5 <- shp[, {
      out <- list(n_peaks = .N, pct_peaks_strand_biased = 100 * mean(strand_p < 0.01))
      for (m in names(mets)) {
        x <- get(m)
        out[[paste0("median_", m)]] <- median(x)
        out[[paste0("p_vs_outside_", m)]] <- if (.BY$category != "Outside CenSat regions" && nrow(ref) >= 5 && .N >= 5)
          suppressWarnings(wilcox.test(x, ref[[m]], exact = FALSE)$p.value) else NA_real_
      }
      out
    }, by = category]
    for (m in names(mets)) s5[, (paste0("q_vs_outside_", m)) := bh(get(paste0("p_vs_outside_", m)))]
    s5[, censat_class := as.character(category)][, category := NULL]
    setcolorder(s5, "censat_class")
    wt(s5[order(-n_peaks)], cfg, "qc05_peak_strand_stacking_by_censat")

    top5 <- head(levels(shp$category), 12)
    lp <- melt(shp[category %in% top5, c("category", names(mets)), with = FALSE], id.vars = "category", variable.name = "metric")
    lp[, metric := factor(wrap(mets[as.character(metric)], 70), levels = wrap(mets, 70))]
    cnt5 <- shp[category %in% top5, .N, by = category]
    lp[, category := factor(category, levels = top5, labels = wrap(paste0(top5, " (n=", cnt5$N[match(top5, cnt5$category)], ")"), 18))]
    p6 <- ggplot(lp, aes(category, value)) +
      geom_violin(data = lp[, if (.N >= 2) .SD, by = .(category, metric)], fill = "#b5d4f4", colour = "#2a78d6", scale = "width", linewidth = 0.3) +
      stat_summary(fun = median, geom = "point", size = 1.6) +
      facet_wrap(~metric, ncol = 1, scales = "free_y") +
      labs(x = NULL, y = NULL, title = "Strand balance and read stacking in peaks, by CenSat class",
           subtitle = sprintf("Peaks with at least %d reads (read 1 of each fragment); black dot = median", minr)) +
      theme_x(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
    save_plot(p6, cfg, "qc05_peak_strand_stacking", max(10, 2 + 0.8 * length(top5)), 10)
    fig_note(cfg, "qc05_peak_strand_stacking",
             sprintf("%d of %d peaks (all IP samples) have at least %d reads and are shown; tests against peaks outside CenSat are in qc05_peak_strand_stacking_by_censat.",
                     nrow(shp), n_all, minr))
  } else fig_note(cfg, "qc05__summary", sprintf("No peak has at least %d reads (qc.shape_min_reads).", minr))
}

## =============================== 5. peak summits relative to CenSat array edges (QC6) ==========
## Mapping artifacts gather at boundaries between unique and repetitive sequence. Observed distances from
## each summit to the nearest edge of its array are compared with random positions in the same arrays.
ap <- rd(file.path(Q, "array_position.tsv"))
if (!is.null(ap) && nrow(ap)) {
  ap <- ap[array_end > array_start + 1]
  ap[, dist := pmin(summit - array_start, array_end - 1 - summit)]
  ap[, array_len := array_end - array_start]
  draws <- qcc$edge_null_draws %||% 20
  set.seed(1)
  nul <- ap[rep(seq_len(.N), each = draws)]
  nul[, pos := array_start + floor(runif(.N) * array_len)]
  nul[, dist := pmin(pos - array_start, array_end - 1 - pos)]
  minp <- qcc$edge_min_peaks %||% 10
  keep <- ap[, .N, by = class][N >= minp, class]
  s6 <- ap[class %in% keep, {
    nd <- nul[class == .BY$class, dist]
    ks <- suppressWarnings(ks.test(dist, nd))
    .(n_peaks = .N, n_arrays = uniqueN(paste(array_start, array_end)), median_array_length = as.numeric(median(array_len)),
      median_distance_observed = as.numeric(median(dist)), median_distance_random = as.numeric(median(nd)),
      pct_within_1kb_observed = 100 * mean(dist < 1000), pct_within_1kb_random = 100 * mean(nd < 1000),
      median_relative_position = as.numeric(median(dist / array_len)), ks_D = unname(ks$statistic), p_value = ks$p.value,
      direction = fifelse(median(dist) < median(nd), "closer to edges than random", "farther from edges than random"))
  }, by = class]
  if (nrow(s6)) {
    s6[, q_value := bh(p_value)]
    s6[, class_name := cls_name(class)]
    setcolorder(s6, c("class", "class_name"))
    wt(s6[order(q_value)], cfg, "qc06_array_edge_distance_by_censat")
    top6 <- head(s6[order(-n_peaks), class], 12)
    lab6 <- s6[class %in% top6, .(class, panel = wrap(sprintf("%s (%d peaks)", class_name, n_peaks), 32),
                                  txt = sprintf("KS D = %.2f, q %s", ks_D, ifelse(q_value < 2.2e-16, fmt_q(q_value), paste("=", fmt_q(q_value)))))]
    ec <- rbind(ap[class %in% top6, .(class, dist, set = "Peak summits")], nul[class %in% top6, .(class, dist, set = "Random positions in the same arrays")])
    ec <- merge(ec, lab6[, .(class, panel)], by = "class")
    lab6[, panel := factor(panel, levels = lab6$panel[match(top6, lab6$class)])]; ec[, panel := factor(panel, levels = levels(lab6$panel))]
    p7 <- ggplot(ec, aes(dist + 1, colour = set)) + stat_ecdf(linewidth = 0.7) +
      geom_text(data = lab6, aes(label = txt), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.6, inherit.aes = FALSE) +
      scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      scale_colour_manual(values = c("Peak summits" = "#e34948", "Random positions in the same arrays" = "grey45"), name = NULL) +
      facet_wrap(~panel, ncol = 4) +
      labs(x = "Distance from summit to the nearest edge of its CenSat array (bp, log scale)", y = "Cumulative fraction",
           title = "Position of peak summits within CenSat arrays") + theme_x(9)
    save_plot(p7, cfg, "qc06_array_edge_distance", 13, 1.8 + 2.6 * ceiling(length(top6) / 4))
    fig_note(cfg, "qc06_array_edge_distance",
             sprintf("%d peak summits in %d classes with at least %d peaks (all samples; a region bound in several samples counts once per sample); %d random positions per peak.",
                     sum(s6$n_peaks), nrow(s6), minp, draws))
  } else fig_note(cfg, "qc06__summary", sprintf("No CenSat class has at least %d peak summits (qc.edge_min_peaks).", minp))
}
