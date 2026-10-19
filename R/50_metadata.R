#!/usr/bin/env Rscript
## Q5: Does sample metadata (virus, cell type, time, cGAS construct, ...) relate to peak load,
## annotation, repeat enrichment, motif content or known-motif (homologous binding site) enrichment?
## Run/batch IDs are never used as model terms or plot keys.
## Tables:  table06_sample_features, table06e_sample_clusters
## Figures: fig16_pca_by_metadata
## (Covariate effects on binding are estimated by the binding models, step 17.)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
cfg <- args_config()
unlink(file.path(cfg$tables, "table06b_confounding.tsv"))   # superseded by table06b_covariate_associations
S <- samples_dt(cfg)
S <- S[peaks_file != ""]
cat_v <- intersect(unlist(cfg$categorical_covariates), names(S))
num_v <- intersect(unlist(cfg$numeric_covariates), names(S))
tv <- cfg$time_var; sv <- cfg$series_var

## readable covariate names (optional config: covariate_labels: {hpi: "Hours post-infection"})
cov_lab <- function(v) {
  l <- cfg$covariate_labels[[v]]
  if (!is.null(l)) return(l)
  if (v == "log10_depth") return("Sequencing depth (log10 mapped reads)")
  gsub("_", " ", v)
}

## ======================= 1. sample x feature matrix ===============================================
feat <- list()
FL <- list()   # feature -> readable label
addlab <- function(f, lab, grp) FL[[length(FL) + 1]] <<- data.table(feature = f, feature_label = lab, feature_group = grp)
ps <- rd(file.path(cfg$work, "peaks", "peak_stats.tsv"))
feat$peaks <- ps[, .(sample_id, peaks__log2_n_peaks = log2(n_peaks + 1), peaks__median_width = median_width)]
addlab("peaks__log2_n_peaks", "Number of peaks (log2)", "Peak load")
addlab("peaks__median_width", "Median peak width (bp)", "Peak load")
addlab("peaks__frip", "Fraction of reads in peaks (FRiP)", "Peak load")

rd_num <- function(f, pattern) {
  if (!nzchar(f) || !file.exists(f)) return(NA_real_)
  l <- grep(pattern, readLines(f), value = TRUE)[1]
  if (is.na(l)) return(NA_real_)
  as.numeric(regmatches(l, regexpr("[0-9.eE+-]+", l)))
}
S[, depth := vapply(flagstat_file, rd_num, numeric(1), pattern = " mapped \\(")]
S[, frip := vapply(frip_file, function(f) {
  if (!nzchar(f) || !file.exists(f)) return(NA_real_)
  l <- grep("^#", readLines(f), value = TRUE, invert = TRUE); l <- l[nzchar(l)]
  if (!length(l)) return(NA_real_)
  suppressWarnings(as.numeric(tail(strsplit(tail(l, 1), "\t")[[1]], 1)))
}, numeric(1))]
feat$qc <- S[, .(sample_id, peaks__frip = frip)]
S[, log10_depth := log10(depth)]

homer <- rbindlist(lapply(seq_len(nrow(S)), function(i) {
  f <- S$homer_file[i]
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  h <- fread(f, sep = "\t", quote = "", select = 8, col.names = "Annotation")
  h[, cat := trimws(sub(" \\(.*$", "", Annotation))]
  h[is.na(cat) | cat == "", cat := "NA"]
  h[, .(pct = 100 * .N / nrow(h)), by = cat][, sample_id := S$sample_id[i]]
}))
if (nrow(homer)) {
  hc <- unique(homer$cat)
  homer[, cat := paste0("homer__pct_", gsub("[^A-Za-z0-9]+", "_", cat))]
  for (x in hc) addlab(paste0("homer__pct_", gsub("[^A-Za-z0-9]+", "_", x)), paste0("% of peaks in ", pretty_homer(x)), "Genomic annotation (HOMER)")
  feat$homer <- dcast(homer, sample_id ~ cat, value.var = "pct", fill = 0)
}

min_obs <- cfg$repeats$min_observed %||% 5
t3 <- rd(file.path(cfg$tables, "table03_repeat_enrichment.tsv"))
if (!is.null(t3)) {
  keep <- t3[, .(mx = max(observed)), by = .(feature_set, label)][mx >= min_obs]
  t3k <- merge(t3, keep[, .(feature_set, label)], by = c("feature_set", "label"))
  t3k[, f := paste0(feature_set, "__oe_", gsub("[^A-Za-z0-9]+", "_", label))]
  u <- unique(t3k[, .(f, feature_set, label)])
  for (i in seq_len(nrow(u))) addlab(u$f[i], paste0("Peak enrichment in ", pretty_label(u$feature_set[i], u$label[i]), " (log2 O/E)"),
                                     if (u$feature_set[i] == "censat") "Centromeric/satellite enrichment" else "RepeatMasker enrichment")
  feat$repeats <- dcast(t3k, sample_id ~ f, value.var = "log2_oe")
}
pv <- rd(file.path(cfg$tables, "table02_motif_prevalence.tsv"))
if (!is.null(pv)) {
  ## families enter the models if detected anywhere OR found de novo in >= q5_min_denovo_samples samples,
  ## keeping the most variable (q5_max_motif_features) so a few samples are not swamped by features
  min_dn <- cfg$motifs$q5_min_denovo_samples %||% 2
  sel <- pv[, .(nd = sum(detected), ndn = sum(discovered_denovo), v = var(pct_peaks)), by = family_id][(nd >= 1 | ndn >= min_dn) & v > 0]
  fk <- head(sel[order(-v)], cfg$motifs$q5_max_motif_features %||% 40)$family_id
  u <- unique(pv[family_id %in% fk, .(family_id, pwm_consensus)])
  for (i in seq_len(nrow(u))) addlab(paste0("motif__", u$family_id[i]), sprintf("Motif family %s (%s): %% of peaks", u$family_id[i], u$pwm_consensus[i]), "Motif family prevalence")
  feat$motifs <- dcast(pv[family_id %in% fk, .(sample_id, f = paste0("motif__", family_id), pct_peaks)], sample_id ~ f, value.var = "pct_peaks")
}
sea <- rd(file.path(cfg$tables, "table04d_sea_known_long.tsv"))
t4 <- rd(file.path(cfg$tables, "table04_known_motif_ranking.tsv"))
if (!is.null(sea) && !is.null(t4)) {
  topk <- head(t4[n_samples_sea > 0], 20)$motif_id
  g <- CJ(sample_id = S$sample_id, motif_id = topk)
  g <- merge(g, sea[, .(sample_id, motif_id, enr_ratio)], by = c("sample_id", "motif_id"), all.x = TRUE)
  g[is.na(enr_ratio), enr_ratio := 1]  # not reported by SEA = not enriched
  g[, f := paste0("known__", gsub("[^A-Za-z0-9]+", "_", t4$tf_name[match(motif_id, t4$motif_id)]), "_", motif_id)]
  u <- unique(g[, .(f, motif_id)])
  for (i in seq_len(nrow(u))) addlab(u$f[i], sprintf("Known motif %s (%s): log2 SEA enrichment", t4$tf_name[match(u$motif_id[i], t4$motif_id)], u$motif_id[i]), "Known-motif enrichment (SEA)")
  feat$known <- dcast(g, sample_id ~ f, value.var = "enr_ratio", fun.aggregate = function(x) log2(x[1]))
}
FM <- Reduce(function(a, b) merge(a, b, by = "sample_id", all = TRUE), feat[!vapply(feat, is.null, logical(1))])
FL <- unique(rbindlist(FL), by = "feature")
meta <- S[, c("sample_id", unique(c(sv, cat_v, num_v, tv))), with = FALSE]
meta[, log10_depth := S$log10_depth]
X <- merge(meta, FM, by = "sample_id")
wt(X, cfg, "table06_sample_features")
fcols <- setdiff(names(FM), "sample_id")
fcols <- fcols[vapply(fcols, function(v) { x <- X[[v]]; sum(!is.na(x)) >= 3 && isTRUE(sd(x, na.rm = TRUE) > 0) }, logical(1))]
flab <- function(f) { l <- FL$feature_label[match(f, FL$feature)]; ifelse(is.na(l), f, l) }
fgrp <- function(f) { g <- FL$feature_group[match(f, FL$feature)]; ifelse(is.na(g), sub("__.*$", "", f), g) }

active_covs <- c(cat_v, num_v)
active_covs <- active_covs[vapply(active_covs, function(v) uniqueN(na.omit(X[[v]])) > 1, logical(1))]

## Earlier versions also wrote per-feature linear models (table06_metadata_effects, Figure 14), linear time-course trends
## (table06c/06d, Figure 15), detection tests (table07), PERMANOVA (table07b) and covariate associations (table06b).
## They were replaced by the binding models (step 17); remove their outputs from earlier runs.
for (t in c("table06_metadata_effects", "table06b_covariate_associations", "table06c_timecourse_trends",
            "table06d_timecourse_series_differences", "table07_motif_detection_vs_metadata", "table07b_permanova"))
  unlink(file.path(cfg$tables, paste0(t, ".tsv")))
clear_figures(cfg, "^fig1[45]_")
for (f in c("fig14", "fig15")) fig_note(cfg, paste0(f, "__summary"), "")

## ======================= 2. PCA and sample clusters (covariates only) ====================================
blocks <- list(`motif prevalence` = grep("^motif__", fcols, value = TRUE),
               `repeat and satellite enrichment` = grep("^(repeatmasker|censat)__", fcols, value = TRUE),
               `genomic annotation` = grep("^homer__", fcols, value = TRUE),
               `known-motif enrichment` = grep("^known__", fcols, value = TRUE))
blocks <- blocks[lengths(blocks) >= 2]
prep <- function(cols) {
  m <- as.matrix(X[, ..cols]); rownames(m) <- X$sample_id
  for (j in seq_len(ncol(m))) m[is.na(m[, j]), j] <- mean(m[, j], na.rm = TRUE)
  m <- m[, apply(m, 2, sd) > 0, drop = FALSE]
  scale(m)
}
allc <- unlist(blocks)
if (length(allc) >= 2 && nrow(X) >= 3 && requireNamespace("patchwork", quietly = TRUE) && length(active_covs)) {
  m <- prep(allc)
  pc <- prcomp(m)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  sc <- data.table(sample_id = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2])
  sc <- merge(sc, X[, c("sample_id", active_covs), with = FALSE], by = "sample_id")
  ## ---- sample clustering on the same scaled features: Ward hierarchical clustering, k by silhouette ----
  D <- dist(m)
  hc <- hclust(D, method = "ward.D2")
  sil_width <- function(cl) {           # mean silhouette width (base R; no extra packages)
    Dm <- as.matrix(D); s <- numeric(length(cl))
    for (i in seq_along(cl)) {
      same <- cl == cl[i]; same[i] <- FALSE
      if (!any(same)) { s[i] <- 0; next }
      a <- mean(Dm[i, same]); b <- min(vapply(setdiff(unique(cl), cl[i]), function(g) mean(Dm[i, cl == g]), 0))
      s[i] <- (b - a) / max(a, b)
    }
    s
  }
  ks <- 2:max(2, min(8, nrow(m) - 1))
  sw <- vapply(ks, function(kk) mean(sil_width(cutree(hc, kk))), 0)
  kbest <- ks[which.max(sw)]
  cl <- cutree(hc, kbest)
  sc[, cluster := factor(paste0("C", cl[sample_id]))]
  sc[, silhouette := sil_width(cl)[match(sample_id, names(cl))]]
  ctab <- merge(sc[, .(sample_id, cluster, silhouette, PC1, PC2)], X[, c("sample_id", active_covs), with = FALSE], by = "sample_id")
  ctab[, n_clusters := kbest][, mean_silhouette := max(sw)]
  wt(ctab[order(cluster, sample_id)], cfg, "table06e_sample_clusters")
  cv_assoc <- vapply(intersect(active_covs, cat_v), function(v) {
    tb <- table(ctab$cluster, ctab[[v]]); kk <- min(dim(tb)) - 1
    if (kk < 1) return(NA_real_)
    sqrt(suppressWarnings(chisq.test(tb, correct = FALSE)$statistic) / (sum(tb) * kk))
  }, 0)
  hull <- sc[, .SD[if (.N >= 3) chull(PC1, PC2) else seq_len(.N)], by = cluster]
  p_cl <- ggplot(sc, aes(PC1, PC2, colour = cluster)) +
    geom_polygon(data = hull, aes(fill = cluster), alpha = 0.12, linewidth = 0.4, show.legend = FALSE) +
    geom_point(size = 2.5) + geom_text(aes(label = short_id(sample_id)), size = 2.2, vjust = -0.9, check_overlap = TRUE, show.legend = FALSE) +
    labs(x = sprintf("PC1 (%s%% of variance)", ve[1]), y = sprintf("PC2 (%s%%)", ve[2]),
         title = sprintf("Sample clusters (Ward clustering, k = %d by silhouette)", kbest), colour = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.2)) + scale_y_continuous(expand = expansion(mult = 0.15)) + theme_x(9)

  plots <- lapply(active_covs, function(v) {
    p <- ggplot(sc, aes(PC1, PC2)) + labs(x = sprintf("PC1 (%s%% of variance)", ve[1]), y = sprintf("PC2 (%s%%)", ve[2]), title = cov_lab(v)) + theme_x(9)
    if (is.numeric(sc[[v]])) p + geom_point(aes(colour = .data[[v]]), size = 2.5) + scale_colour_gradient(low = "#b5d4f4", high = "#0c447c", name = NULL)
    else p + geom_point(aes(colour = factor(.data[[v]]), shape = factor(.data[[v]])), size = 2.5) + labs(colour = NULL, shape = NULL)
  })
  pp <- patchwork::wrap_plots(c(list(p_cl), plots), ncol = 2)
  save_plot(pp, cfg, "fig16_pca_by_metadata", 12, 4.5 * ceiling((length(plots) + 1) / 2))
  fig_note(cfg, "fig16_pca_by_metadata", sprintf("PCA of %d scaled sample features (%s). First panel: %d clusters (mean silhouette %.2f; k tested 2-%d). Agreement of clusters with covariates (Cramer's V): %s.",
           ncol(m), paste(names(blocks), collapse = ", "), kbest, max(sw), max(ks),
           if (length(cv_assoc)) paste(sprintf("%s %.2f", vapply(names(cv_assoc), cov_lab, ""), cv_assoc), collapse = ", ") else "none"))
}
