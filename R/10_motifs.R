#!/usr/bin/env Rscript
## Q1: How often does each motif occur per sample, and does it look like a repeat?
## Tables:  table01_master_motifs, table01b_motif_family_members, table02_motif_prevalence (long),
##          table02_motif_prevalence_matrix, table02_repeat_likeness, table02c_dfam_matches
## Figures: fig01_motif_prevalence_heatmap, fig02_repeat_likeness_scatter, fig03_motif_spacing
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
suppressPackageStartupMessages(library(pheatmap))
cfg <- args_config()
mc <- cfg$motifs
S <- samples_dt(cfg)

fam <- rd(file.path(cfg$work, "motifs", "families.tsv"))
feat <- rd(file.path(cfg$work, "motifs", "family_rep_features.tsv"))
xm <- rd(file.path(cfg$tables, "xstreme_motifs.tsv"))
info <- rd(file.path(cfg$tables, "xstreme_run_info.tsv"))
if (is.null(fam) || is.null(feat)) stop("run steps 02-04 first")

## ---------------- master motif table (Table 1) -----------------------------------------------------
kidx <- rd(file.path(cfg$work, "known", "known_motif_index.tsv"))
ttk <- read_tomtom(file.path(cfg$work, "known", "tomtom_discovered_vs_known.tsv"))
best_known <- NULL
if (!is.null(ttk) && nrow(ttk)) {
  best_known <- ttk[e_value <= (mc$known_evalue %||% 0.05)][order(e_value)][, .SD[1], by = query]
  if (!is.null(kidx)) best_known <- merge(best_known, kidx[, .(target = motif_id, target_name = alt_id)], by = "target", all.x = TRUE)
  else best_known[, target_name := target]
  best_known <- best_known[, .(motif_uid = query, best_known_match = paste0(target_name, " (", target, ")"),
                               best_known_evalue = e_value)]
}
master <- merge(fam[, .(motif_uid, family_id, is_representative, family_n_samples)],
                xm[is_denovo == TRUE], by = "motif_uid", all.x = TRUE)
master[, sample_id := sub("__.*$", "", motif_uid)]
if (!is.null(best_known)) master <- merge(master, best_known, by = "motif_uid", all.x = TRUE)
master <- merge(master, S[, .SD, .SDcols = intersect(names(S), c("sample_id", cfg$batch_var, unlist(cfg$categorical_covariates), unlist(cfg$numeric_covariates)))],
                by = "sample_id", all.x = TRUE)
setcolorder(master, c("family_id", "sample_id", "motif_uid"))
setorder(master, family_id, sample_id)
wt(master, cfg, "table01_master_motifs")
if (!is.null(info) && any(info$known_db_used == FALSE & info$status == "ok"))
  msg("NOTE: ", sum(info$known_db_used == FALSE & info$status == "ok"),
      " XSTREME runs had no known-motif DB; their sim_motif columns are empty (Tomtom here fills the gap)")

## ---------------- prevalence from the uniform FIMO rescans ------------------------------------------
sm <- rbindlist(lapply(list.files(file.path(cfg$work, "scan"), "\\.(peaks|background)\\.summary\\.tsv$", full.names = TRUE), rd))
if (!nrow(sm)) stop("no scan summaries; run step 05")
pk <- sm[set == "peaks", .(sample_id, family_id, n_seqs, n_hit = n_seqs_hit, n_sites, max_sites = max_sites_per_seq)]
bg <- sm[set == "background", .(sample_id, family_id, bg_n = n_seqs, bg_hit = n_seqs_hit)]
prev <- merge(pk, bg, by = c("sample_id", "family_id"), all.x = TRUE)
prev[, pct_peaks := 100 * n_hit / pmax(n_seqs, 1)]
prev[, pct_background := 100 * bg_hit / pmax(bg_n, 1)]
## fold with a half-count correction, so a handful of hits against 0 background hits cannot give
## absurd ratios (the previous tiny pseudocount turned e.g. 2 vs 0 hits into ~170-fold)
prev[, enrichment := ((n_hit + 0.5) / (n_seqs + 1)) / ((bg_hit + 0.5) / (bg_n + 1))]
prev[, mean_sites_per_hit_peak := fifelse(n_hit > 0, n_sites / n_hit, NA_real_)]
prev[, p_fisher := mapply(function(a, n, b, m) {
  if (n == 0 || m == 0) return(NA_real_)
  fisher.test(matrix(c(a, n - a, b, m - b), 2), alternative = "greater")$p.value
}, n_hit, n_seqs, bg_hit, bg_n)]
prev[, q_fisher := bh(p_fisher)]
dq <- mc$detect_qvalue %||% 0.05; de <- mc$detect_min_enrichment %||% 1.5
prev[, detected_vs_genomic := !is.na(q_fisher) & q_fisher < dq & enrichment >= de]

## SEA of each family vs shuffled copies of the sample's own peaks (same kind of control XSTREME used)
fsea <- rbindlist(lapply(list.files(file.path(cfg$work, "scan"), "\\.family_sea\\.tsv$", full.names = TRUE), function(f) {
  l <- readLines(f); l <- l[!grepl("^#", l) & nzchar(l)]
  if (length(l) < 2) return(NULL)
  d <- fread(text = paste(l, collapse = "\n"))
  setnames(d, tolower(gsub("%", "_pct", gsub("[^A-Za-z0-9%]+", "_", names(d)))))
  d[, .(sample_id = sub("\\.family_sea\\.tsv$", "", basename(f)), family_id = id,
        sea_enr_ratio = as.numeric(enr_ratio), sea_qvalue = as.numeric(qvalue))]
}), fill = TRUE)
if (nrow(fsea)) prev <- merge(prev, fsea, by = c("sample_id", "family_id"), all.x = TRUE) else
  prev[, `:=`(sea_enr_ratio = NA_real_, sea_qvalue = NA_real_)]
prev[, detected_vs_shuffled := !is.na(sea_qvalue) & sea_qvalue < dq & sea_enr_ratio >= de]
against <- mc$detect_against %||% "genomic"
prev[, detected := switch(against, genomic = detected_vs_genomic, shuffled = detected_vs_shuffled,
                          both = detected_vs_genomic & detected_vs_shuffled, either = detected_vs_genomic | detected_vs_shuffled,
                          stop("motifs.detect_against must be genomic, shuffled, both or either"))]

## diagnosis printed to the log: why are families (not) detected?
q50 <- function(x) sprintf("%.3g", median(x, na.rm = TRUE))
msg(sprintf("detection diagnosis over %d family x sample tests: median %% peaks with motif %s, median %% background %s, median fold %s; ",
            nrow(prev), q50(prev$pct_peaks), q50(prev$pct_background), q50(prev$enrichment)),
    sprintf("vs genomic: %d with q < %s, %d with fold >= %s, %d both; vs shuffled (SEA): %d detected; using '%s'.",
            sum(prev$q_fisher < dq, na.rm = TRUE), dq, sum(prev$enrichment >= de, na.rm = TRUE), de,
            sum(prev$detected_vs_genomic), sum(prev$detected_vs_shuffled), against))
if (median(prev$pct_background, na.rm = TRUE) > 50)
  msg("NOTE: most background windows contain a site, so prevalence saturates; use seq_mode: summit and/or a smaller fimo_thresh.")
denovo <- unique(fam[, .(sample_id, family_id, discovered_denovo = TRUE)])
prev <- merge(prev, denovo, by = c("sample_id", "family_id"), all.x = TRUE)
prev[is.na(discovered_denovo), discovered_denovo := FALSE]
prev <- merge(prev, feat[, .(family_id, pwm_consensus)], by = "family_id", all.x = TRUE)
setcolorder(prev, c("family_id", "pwm_consensus", "sample_id"))
wt(prev, cfg, "table02_motif_prevalence")
wt(dcast(prev, family_id + pwm_consensus ~ sample_id, value.var = "pct_peaks"), cfg, "table02_motif_prevalence_matrix")

## ---------------- repeat-likeness (Table 2) ----------------------------------------------------------
sp <- rbindlist(lapply(list.files(file.path(cfg$work, "scan"), "\\.peaks\\.spacing\\.tsv$", full.names = TRUE), rd))
per <- if (nrow(sp)) sp[spacing >= 2, .(n = sum(n)), by = .(family_id, spacing)] else data.table(family_id = character(), spacing = integer(), n = numeric())
spacing_stats <- per[, {
  mode <- spacing[which.max(n)]
  tot <- sum(n)
  pf <- if (mode >= 3) sum(n[(spacing %% mode) %in% c(0, 1, mode - 1)]) / tot else NA_real_
  .(spacing_mode = mode, spacing_mode_frac = max(n) / tot, periodic_frac = pf, n_spacings = tot)
}, by = family_id]

sl <- rbindlist(lapply(list.files(file.path(cfg$work, "scan"), "\\.site_labels\\.tsv$", full.names = TRUE), rd))
tot_sites <- pk[, .(total_sites = sum(n_sites)), by = family_id]
rep_stats <- data.table(family_id = feat$family_id)
if (nrow(sl)) {
  anyl <- sl[label == "ANY", .(n = sum(n_sites)), by = .(family_id, feature_set)]
  anyl <- dcast(anyl, family_id ~ feature_set, value.var = "n", fill = 0)
  top <- sl[label != "ANY", .(n = sum(n_sites)), by = .(family_id, feature_set, label)][order(-n)][, .SD[1], by = family_id]
  rep_stats <- merge(rep_stats, anyl, by = "family_id", all.x = TRUE)
  rep_stats <- merge(rep_stats, top[, .(family_id, top_repeat_label = paste0(feature_set, ":", label), top_label_sites = n)], by = "family_id", all.x = TRUE)
}
rep_stats <- merge(rep_stats, tot_sites, by = "family_id", all.x = TRUE)
for (v in c("repeatmasker", "censat", "top_label_sites")) if (!v %in% names(rep_stats)) rep_stats[, (v) := NA_real_]
if (!"top_repeat_label" %in% names(rep_stats)) rep_stats[, top_repeat_label := NA_character_]
rep_stats[, `:=`(pct_sites_in_repeatmasker = 100 * repeatmasker / pmax(total_sites, 1),
                 pct_sites_in_censat = 100 * censat / pmax(total_sites, 1),
                 pct_sites_top_label = 100 * top_label_sites / pmax(total_sites, 1))]

rc <- NULL
rcf <- file.path(cfg$work, "motifs", "reps_vs_repeat_consensus.tsv")
if (file.exists(rcf) && file.info(rcf)$size > 0) {
  rcraw <- fread(rcf)
  setnames(rcraw, tolower(gsub("[^A-Za-z0-9]+", "_", names(rcraw))))

  ## ---- table02c: every Dfam consensus match, with a Benjamini-Hochberg correction ----------------
  ## step 05 runs `fimo --text`, and MEME's documentation states that text mode "will also disable the
  ## calculation and printing of q-values", so FIMO gives no q-value here. The correction below is
  ## applied afterwards, to the matches FIMO actually reported.
  ##
  ## IMPORTANT, and stated again in the table's own columns: FIMO only reports matches passing
  ## motifs.fimo_thresh, so this corrects over the SURVIVORS of that threshold, not over every position
  ## tested. It is therefore anti-conservative relative to a FIMO q-value from a non-text run. Treat it
  ## as a ranking aid among the reported hits, not as a calibrated FDR.
  ##
  ## Two units of correction are given because positional matches within one consensus sequence are not
  ## independent -- a satellite consensus contains the motif many times by construction:
  ##   q_bh_pair  BH over the best match of each family x consensus-sequence pair  (the default; the
  ##              unit of inference is "does this family match this consensus at all")
  ##   q_bh_match BH over every reported positional match
  dm <- copy(rcraw)
  setnames(dm, "motif_id", "family_id")
  if (!"q_value" %in% names(dm)) dm[, q_value := NA_real_]   # absent in --text output
  dm[, q_bh_match := p.adjust(p_value, method = "BH")]
  ## best match per family x sequence, then BH over those pairs
  setorder(dm, family_id, sequence_name, p_value)
  cnt <- dm[, .(n_matches_in_sequence = .N), by = .(family_id, sequence_name)]
  pairs <- merge(dm[, .SD[1], by = .(family_id, sequence_name)], cnt,
                 by = c("family_id", "sequence_name"))       # join, not row-order alignment
  pairs[, q_bh_pair := p.adjust(p_value, method = "BH")]
  dfam <- pairs[, .(family_id, dfam_sequence = sequence_name, best_p = p_value,
                    q_bh_pair, q_bh_match, fimo_q = q_value, n_matches_in_sequence,
                    start, stop, strand, score,
                    n_pairs_corrected = nrow(pairs), n_matches_corrected = nrow(dm),
                    fimo_thresh = mc$fimo_thresh %||% 1e-4,
                    correction = "BH over FIMO-reported matches only; anti-conservative (fimo --text gives no q-value)")]
  setorder(dfam, family_id, best_p)
  dfam[, rank_in_family := seq_len(.N), by = family_id]
  wt(dfam, cfg, "table02c_dfam_matches")
  msg("Dfam matches: ", nrow(dm), " reported matches over ", nrow(pairs),
      " family x consensus pairs; BH applied to both units")

  ## ---- the single best hit per family, as before, now carrying its corrected value ---------------
  rc <- dfam[order(best_p)][, .(repeat_consensus_hit = dfam_sequence[1], repeat_consensus_p = best_p[1],
                                repeat_consensus_q = q_bh_pair[1],
                                n_repeat_consensus_hit = uniqueN(dfam_sequence)), by = family_id]
}

pooled <- prev[, .(samples_discovered_denovo = sum(discovered_denovo), samples_detected = sum(detected),
                   n_samples = .N, median_pct_peaks = median(pct_peaks),
                   mean_sites_per_hit_peak = sum(n_sites) / pmax(sum(n_hit), 1),
                   max_sites_in_one_peak = max(max_sites)), by = family_id]
t2 <- Reduce(function(a, b) merge(a, b, by = "family_id", all.x = TRUE),
             list(feat, pooled, spacing_stats, rep_stats[, .(family_id, total_sites, pct_sites_in_repeatmasker,
                  pct_sites_in_censat, top_repeat_label, pct_sites_top_label)]))
if (!is.null(rc)) t2 <- merge(t2, rc, by = "family_id", all.x = TRUE)
t2[, pct_sites_in_repeats := pmax(pct_sites_in_repeatmasker, pct_sites_in_censat, na.rm = TRUE)]
thr_sites <- mc$repeat_min_sites_per_peak %||% 3
thr_per <- mc$repeat_min_periodic_frac %||% 0.3
thr_rep <- mc$repeat_min_pct_sites_in_repeats %||% 50
## periodicity only counts with enough site pairs to estimate it
t2[, tandem_evidence := (mean_sites_per_hit_peak >= thr_sites) |
       (!is.na(periodic_frac) & periodic_frac >= thr_per & n_spacings >= 30 & mean_sites_per_hit_peak >= 1.5) |
       (self_period_corr >= 0.8 & width >= 8)]
t2[, flag := fifelse(!is.na(pct_sites_in_repeats) & pct_sites_in_repeats >= thr_rep & tandem_evidence, "repeat (tandem)",
             fifelse(!is.na(pct_sites_in_repeats) & pct_sites_in_repeats >= thr_rep, "repeat-embedded",
             fifelse(tandem_evidence | trimer_diversity < 0.5, "low-complexity / periodic", "TF-like")))]
setorder(t2, -samples_detected, -samples_discovered_denovo, -median_pct_peaks)
wt(t2, cfg, "table02_repeat_likeness")

## ---------------- figures ---------------------------------------------------------------------------
## Fig 1: prevalence heatmap (no row-annotation key; sample metadata only, no run ID).
## Rows: the 5 families most enriched over the genomic background in each sample (fold enrichment, then % of peaks),
## combined over samples, so common repeat motifs do not take every slot.
top_n1 <- 5
top_f <- unique(prev[n_hit > 0][order(-enrichment, -pct_peaks), head(family_id, top_n1), by = sample_id]$V1)
mat <- dcast(prev[family_id %in% top_f], family_id ~ sample_id, value.var = "pct_peaks")
m <- as.matrix(mat[, -1]); rownames(m) <- paste(mat$family_id, feat$pwm_consensus[match(mat$family_id, feat$family_id)])
if (nrow(m) >= 1 && ncol(m) >= 1) {
  ann_c <- sample_annotation(cfg, S, colnames(m))
  save_pheatmap(function() pheatmap(m, cluster_rows = nrow(m) > 2, cluster_cols = ncol(m) > 2,
                                    annotation_col = ann_c, fontsize = 8,
                                    color = colorRampPalette(c("white", "#2a78d6", "#0c447c"))(50),
                                    main = "Motif family prevalence (% of peaks with >= 1 site)"),
                cfg, "fig01_motif_prevalence_heatmap", width = max(8, 3 + 0.3 * ncol(m)), height = max(5, 3 + 0.22 * nrow(m)))
  fig_note(cfg, "fig01_motif_prevalence_heatmap",
           sprintf("Showing %d of %d motif families: the %d most enriched over the genomic background (fold enrichment, then %% of peaks) in each of %d samples, combined.",
                   nrow(m), nrow(t2), top_n1, ncol(m)))
}

## Fig 2: repeat-likeness; label only outliers (beyond 1.5 x IQR on either axis, log scale for y)
iqr_out <- function(x) { q <- quantile(x, c(.25, .75), na.rm = TRUE); !is.na(x) & (x > q[2] + 1.5 * diff(q) | x < q[1] - 1.5 * diff(q)) }
t2p <- copy(t2[!is.na(mean_sites_per_hit_peak) & mean_sites_per_hit_peak > 0])   # 0 sites cannot go on a log axis
t2p[, is_outlier := iqr_out(median_pct_peaks) | iqr_out(log10(mean_sites_per_hit_peak))]
p2 <- ggplot(t2p, aes(median_pct_peaks, mean_sites_per_hit_peak, colour = pct_sites_in_repeats, shape = flag)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(data = t2p[is_outlier == TRUE], aes(label = family_id), size = 3, vjust = -1, show.legend = FALSE, colour = "grey20") +
  scale_colour_gradient(low = "#2a78d6", high = "#eb6834", limits = c(0, 100), na.value = "grey70",
                        name = "% sites in repeats") +
  scale_y_log10() +
  labs(x = "Median % of peaks with the motif", y = "Mean sites per occupied peak (log scale)",
       title = "Repeat-likeness of motif families", shape = NULL) + theme_x()
save_plot(p2, cfg, "fig02_repeat_likeness_scatter", 8, 6)
fig_note(cfg, "fig02_repeat_likeness_scatter",
         sprintf("%d motif families plotted (%d with no sites in any peak omitted); %d outliers labelled.",
                 nrow(t2p), nrow(t2) - nrow(t2p), sum(t2p$is_outlier)))

## Fig 3: spacing histograms for the 20 families with the strongest periodicity
if (nrow(per)) {
  rk <- t2[!is.na(spacing_mode_frac) & n_spacings >= 30][order(-spacing_mode_frac, -n_spacings)]
  if (!nrow(rk)) rk <- t2[!is.na(n_spacings)][order(-n_spacings)]
  showf <- head(rk$family_id, 20)
  ps <- per[family_id %in% showf]
  ps[, family := factor(paste(family_id, feat$pwm_consensus[match(family_id, feat$family_id)]),
                        levels = paste(showf, feat$pwm_consensus[match(showf, feat$family_id)]))]
  p3 <- ggplot(ps, aes(spacing, n)) + geom_col(width = 1, fill = "#2a78d6") +
    facet_wrap(~family, scales = "free_y", ncol = 4) +
    labs(x = "Distance between adjacent same-strand sites (bp)", y = "Site pairs",
         title = "Motif site spacing (fixed periods indicate tandem repeats)") + theme_x(9)
  save_plot(p3, cfg, "fig03_motif_spacing", 13, 1.5 + 2.2 * ceiling(length(showf) / 4))
  fig_note(cfg, "fig03_motif_spacing",
           sprintf("Top %d of %d families ranked by periodicity strength (share of site pairs at the most common spacing; families with >= 30 site pairs).",
                   length(showf), nrow(t2)))
}
