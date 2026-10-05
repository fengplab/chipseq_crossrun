#!/usr/bin/env Rscript
## Q5: Does sample metadata (virus, cell type, time, cGAS construct, ...) relate to peak load,
## annotation, repeat enrichment, motif content or known-motif (homologous binding site) enrichment?
## Run/batch IDs are never used as model terms or plot keys.
## Tables:  table06_sample_features, table06_metadata_effects, table06b_covariate_associations,
##          table06c_timecourse_trends, table06d_timecourse_series_differences,
##          table07_motif_detection_vs_metadata, table07b_permanova
## Figures: fig14_effects_heatmap, fig15_timecourse_<series>, fig16_pca_by_metadata
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
cfg <- args_config()
unlink(file.path(cfg$tables, "table06b_confounding.tsv"))   # superseded by table06b_covariate_associations
S <- samples_dt(cfg)
S <- S[peaks_file != ""]
cat_v <- intersect(unlist(cfg$categorical_covariates), names(S))
num_v <- intersect(unlist(cfg$numeric_covariates), names(S))
tv <- cfg$time_var; sv <- cfg$series_var
tc_cfg <- cfg$timecourse %||% list()
tc_alpha <- tc_cfg$alpha %||% 0.05
tc_use <- tc_cfg$significance %||% "q"   # q (BH within series) or p (nominal)

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

## ======================= 2. how strongly covariates are tied to each other ======================
## (replaces the old run-ID confounding check). Categorical-categorical: Cramer's V;
## categorical-numeric: eta^2 (R^2 of numeric ~ categorical); numeric-numeric: |Spearman rho|.
active_covs <- c(cat_v, num_v)
active_covs <- active_covs[vapply(active_covs, function(v) uniqueN(na.omit(X[[v]])) > 1, logical(1))]
assoc <- rbindlist(if (length(active_covs) >= 2) combn(active_covs, 2, simplify = FALSE, FUN = function(p) {
  d <- X[!is.na(get(p[1])) & !is.na(get(p[2])), p, with = FALSE]
  a <- d[[1]]; b <- d[[2]]
  if (p[1] %in% cat_v && p[2] %in% cat_v) {
    tb <- table(a, b); chi <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
    k <- min(dim(tb)) - 1
    v <- if (k > 0) sqrt(unname(chi) / (sum(tb) * k)) else NA_real_; m <- "Cramer's V"
  } else if (xor(p[1] %in% cat_v, p[2] %in% cat_v)) {
    cc <- if (p[1] %in% cat_v) a else b; nn <- if (p[1] %in% cat_v) b else a
    v <- summary(lm(nn ~ factor(cc)))$r.squared; m <- "eta squared"
  } else { v <- abs(suppressWarnings(cor(a, b, method = "spearman"))); m <- "|Spearman rho|" }
  data.table(covariate_a = p[1], covariate_b = p[2], measure = m, association = v, n = nrow(d),
             note = fifelse(!is.na(v) & v >= 0.99, "aliased: effects cannot be separated",
                    fifelse(!is.na(v) & v >= 0.7, "strongly associated: interpret effects with caution", "")))
}) else list())
if (nrow(assoc)) wt(assoc[order(-association)], cfg, "table06b_covariate_associations")

## ======================= 3. per-feature linear models (no run/batch term) ========================
term_labeler <- function(d, rhs) {
  out <- list()
  for (v in rhs) {
    if (is.factor(d[[v]])) { lv <- levels(d[[v]])
      for (l in lv[-1]) out[[length(out) + 1]] <- data.table(term = paste0(v, l), covariate = v,
                                                              term_label = sprintf("%s: %s vs %s", cov_lab(v), l, lv[1]))
    } else out[[length(out) + 1]] <- data.table(term = v, covariate = v, term_label = sprintf("%s (per 1-unit increase)", cov_lab(v)))
  }
  rbindlist(out)
}
fit_one <- function(y, covs, extra = NULL) {
  rhs <- c(covs, extra)
  d <- X[, c(y, rhs), with = FALSE]
  d <- d[complete.cases(d)]
  for (v in intersect(rhs, cat_v)) d[[v]] <- factor(d[[v]])
  rhs <- rhs[vapply(rhs, function(v) uniqueN(d[[v]]) > 1, logical(1))]
  if (!length(rhs) || nrow(d) < length(rhs) + 3) return(NULL)
  m <- lm(as.formula(paste(sprintf("`%s`", y), "~", paste(sprintf("`%s`", rhs), collapse = " + "))), data = d)
  co <- summary(m)$coefficients
  res <- data.table(term = gsub("`", "", rownames(co)), estimate = co[, 1], std_error = co[, 2], p_value = co[, 4])
  tl <- term_labeler(d, rhs)
  res <- merge(tl, res, by = "term", all.x = TRUE)   # aliased terms keep NA estimates
  res[, `:=`(n = nrow(d), model = "linear model (OLS)")][]
}
eff <- rbindlist(lapply(fcols, function(y) {
  extra <- if (startsWith(y, "peaks__") && sum(!is.na(X$log10_depth)) >= 3) "log10_depth" else NULL
  r <- fit_one(y, active_covs, extra); mode <- "joint"
  if (is.null(r)) { r <- rbindlist(lapply(active_covs, function(v) fit_one(y, v, extra))); mode <- "one covariate at a time" }
  if (!nrow(r)) return(NULL)
  r[, `:=`(feature = y, feature_label = flab(y), feature_group = fgrp(y), fit = mode)]
}), fill = TRUE)
if (nrow(eff)) {
  eff[, `:=`(ci_low = estimate - 1.96 * std_error, ci_high = estimate + 1.96 * std_error)]
  eff[, note := fifelse(is.na(estimate), "not estimable (aliased with another covariate)", "")]
  eff[, q_value := bh(p_value), by = feature_group]
  setcolorder(eff, c("feature_group", "feature", "feature_label", "covariate", "term", "term_label", "estimate", "std_error",
                     "ci_low", "ci_high", "p_value", "q_value", "model", "fit", "n", "note"))
  setorder(eff, q_value, p_value, na.last = TRUE)
  wt(eff, cfg, "table06_metadata_effects")

  hm <- eff[covariate != "log10_depth" & !is.na(p_value)]
  topf <- head(unique(hm[order(q_value, p_value)]$feature_label), 30)
  hm <- hm[feature_label %in% topf]
  hm[, score := sign(estimate) * pmin(-log10(q_value), 10)]
  p14 <- ggplot(hm, aes(term_label, factor(feature_label, levels = rev(topf)), fill = score)) + geom_tile(colour = "white") +
    geom_text(data = hm[q_value < 0.05], label = "*", vjust = 0.75) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", name = "sign x -log10 q") +
    scale_x_discrete(labels = function(x) stringr_wrap(x, 28)) +
    labs(x = NULL, y = NULL, title = "Metadata effects on sample features", subtitle = "* q < 0.05 within feature group") +
    theme_x(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right")
  save_plot(p14, cfg, "fig14_effects_heatmap", max(10, 6 + 0.7 * uniqueN(hm$term_label)), 3 + 0.3 * length(topf))
  fig_note(cfg, "fig14_effects_heatmap", sprintf("Top %d of %d features, ranked by their smallest q-value across terms.", length(topf), uniqueN(eff$feature)))
}

## ======================= 4. time courses: one figure per series ==================================
clear_figures(cfg, "^fig15_timecourse")
if (!is.null(sv) && !is.null(tv) && all(c(sv, tv) %in% names(X))) {
  tc <- X[!is.na(get(sv)) & get(sv) != "" & !is.na(get(tv))]
  ok_series <- tc[, uniqueN(get(tv)), by = c(sv)][V1 >= 3][[sv]]
  skipped <- setdiff(unique(tc[[sv]]), ok_series)
  tc <- tc[get(sv) %in% ok_series]
  tr <- rbindlist(lapply(fcols, function(y) rbindlist(lapply(ok_series, function(s) {
    d <- tc[get(sv) == s & !is.na(get(y))]
    if (nrow(d) < 3 || uniqueN(d[[y]]) < 2) return(NULL)
    fit <- summary(lm(d[[y]] ~ d[[tv]]))$coefficients
    rho <- suppressWarnings(cor(d[[tv]], d[[y]], method = "spearman"))
    data.table(feature = y, feature_label = flab(y), feature_group = fgrp(y), series = s, n_timepoints = nrow(d),
               slope_per_unit_time = fit[2, 1], p_value = if (nrow(fit) > 1 && ncol(fit) >= 4) fit[2, 4] else NA_real_,
               spearman_rho = rho)
  }))))
  sig_notes <- character()
  if (nrow(tr)) {
    tr[, q_value := bh(p_value), by = series]
    tr[, significant := !is.na(p_value) & (if (tc_use == "p") p_value else q_value) < tc_alpha]
    wt(tr[order(series, p_value)], cfg, "table06c_timecourse_trends")
    for (s in ok_series) {
      sig <- tr[series == s & significant == TRUE][order(p_value)]
      if (!nrow(sig)) { sig_notes <- c(sig_notes, sprintf("%s: no significant relations (%d time points; %s < %s)", s, tc[get(sv) == s, .N], tc_use, tc_alpha)); next }
      show <- head(sig$feature, 20)
      pl <- melt(tc[get(sv) == s, c(tv, show), with = FALSE], id.vars = tv, variable.name = "feature")
      setnames(pl, tv, "time")
      pl[, panel := factor(stringr_wrap(flab(as.character(feature)), 38), levels = stringr_wrap(flab(show), 38))]
      lab <- sig[match(show, feature), sprintf("%s = %.2g", tc_use, if (tc_use == "p") p_value else q_value)]
      ann <- data.table(panel = factor(stringr_wrap(flab(show), 38), levels = stringr_wrap(flab(show), 38)), lab = lab)
      p15 <- ggplot(pl, aes(time, value)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey60", linewidth = 0.5) +
        geom_line(colour = "#2a78d6", linewidth = 0.8) + geom_point(colour = "#2a78d6", size = 2) +
        geom_text(data = ann, aes(label = lab), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
        facet_wrap(~panel, scales = "free_y", ncol = 4) +
        labs(x = cov_lab(tv), y = "Feature value (units in panel title)", title = sprintf("Time course: %s", s)) + theme_x(8)
      fn <- paste0("fig15_timecourse_", gsub("[^A-Za-z0-9]+", "_", s))
      save_plot(p15, cfg, fn, 13, 1.5 + 2.4 * ceiling(length(show) / 4))
      fig_note(cfg, fn, sprintf("%d of %d features change significantly with %s in this series (%s < %s); %s shown.",
                                nrow(sig), tr[series == s, .N], cov_lab(tv), tc_use, tc_alpha,
                                if (nrow(sig) > 20) "top 20 by p-value" else "all"))
    }
  }
  if (length(skipped)) sig_notes <- c(sig_notes, sprintf("%s: not tested (fewer than 3 time points)", skipped))
  fig_note(cfg, "fig15_timecourse__summary", paste(sig_notes, collapse = "; "))
  if (length(ok_series) >= 2) {
    dif <- rbindlist(lapply(fcols, function(y) {
      d <- tc[!is.na(get(y))]
      d[, series_f := factor(get(sv))]
      if (uniqueN(d$series_f) < 2) return(NULL)
      m0 <- lm(as.formula(sprintf("`%s` ~ `%s` + series_f", y, tv)), d)
      m1 <- lm(as.formula(sprintf("`%s` ~ `%s` * series_f", y, tv)), d)
      a <- anova(m0, m1)
      data.table(feature = y, feature_label = flab(y), feature_group = fgrp(y), n = nrow(d), F = a$F[2], p_slope_differs = a$`Pr(>F)`[2])
    }))
    if (nrow(dif)) { dif[, q_value := bh(p_slope_differs), by = feature_group]; wt(dif[order(p_slope_differs)], cfg, "table06d_timecourse_series_differences") }
  }
}

## ======================= 5. motif / known-motif detection vs metadata ============================
det <- list()
if (!is.null(pv)) det$`motif family` <- pv[, .(sample_id, id = paste(family_id, pwm_consensus), detected)]
if (!is.null(sea) && !is.null(t4)) {
  sq <- cfg$motifs$sea_qvalue %||% 0.05
  g <- CJ(sample_id = S$sample_id, motif_id = head(t4[n_samples_sea > 0], 50)$motif_id)
  g <- merge(g, sea[, .(sample_id, motif_id, qvalue)], by = c("sample_id", "motif_id"), all.x = TRUE)
  det$`known motif` <- g[, .(sample_id, id = paste(t4$tf_name[match(motif_id, t4$motif_id)], motif_id), detected = !is.na(qvalue) & qvalue <= sq)]
}
t7 <- rbindlist(lapply(names(det), function(k) {
  d <- merge(det[[k]], X[, c("sample_id", active_covs), with = FALSE], by = "sample_id")
  d <- d[, if (any(detected) && !all(detected)) .SD, by = id]
  if (!nrow(d)) return(NULL)
  rbindlist(lapply(active_covs, function(v) {
    ## test only motifs whose detection varies among the samples that HAVE this covariate: variation that
    ## comes only from samples missing the covariate leaves one comparison group empty
    dv <- d[!is.na(get(v))][, if (any(detected) && !all(detected)) .SD, by = id]
    if (!nrow(dv)) return(NULL)
    d <- dv
    if (v %in% cat_v) {
      rbindlist(lapply(unique(na.omit(d[[v]])), function(lv) d[!is.na(get(v)), {
        inl <- get(v) == lv
        ft <- fisher.test(matrix(c(sum(detected & inl), sum(!detected & inl), sum(detected & !inl), sum(!detected & !inl)), 2))
        .(covariate = v, level = as.character(lv), detected_in_level = sprintf("%d/%d", sum(detected & inl), sum(inl)),
          detected_elsewhere = sprintf("%d/%d", sum(detected & !inl), sum(!inl)),
          odds_ratio = unname(ft$estimate), p_value = ft$p.value, test = "Fisher (level vs rest)")
      }, by = id]))
    } else {
      d[!is.na(get(v)), {
        x1 <- get(v)[detected]; x0 <- get(v)[!detected]
        w <- suppressWarnings(wilcox.test(x1, x0, exact = FALSE))
        .(covariate = v, level = "numeric", detected_in_level = sprintf("median %s = %.3g (detected)", v, median(x1)),
          detected_elsewhere = sprintf("median %s = %.3g (not)", v, median(x0)), odds_ratio = NA_real_,
          p_value = w$p.value, test = "Wilcoxon")
      }, by = id]
    }
  }))[, feature_type := k]
}), fill = TRUE)
if (!nrow(t7)) {
  unlink(file.path(cfg$tables, "table07_motif_detection_vs_metadata.tsv"))
  msg("NOTE: table07 not written: no motif family or known motif is detected in some samples but not others ",
      "(detection is all-or-none for every motif), so there is nothing to associate with metadata.")
} else {
  t7[, q_value := bh(p_value), by = feature_type]
  setcolorder(t7, c("feature_type", "id", "covariate", "level"))
  wt(t7[order(q_value, p_value)], cfg, "table07_motif_detection_vs_metadata")
}

## ======================= 6. PERMANOVA + PCA (covariates only) ====================================
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
if (requireNamespace("vegan", quietly = TRUE) && length(active_covs)) {
  perm <- rbindlist(lapply(names(blocks), function(b) {
    m <- prep(blocks[[b]])
    d <- as.data.frame(X[, active_covs, with = FALSE]); rownames(d) <- X$sample_id
    ok <- complete.cases(d)
    if (sum(ok) < 5 || ncol(m) < 2) return(NULL)
    d <- d[ok, , drop = FALSE]; mm <- m[ok, , drop = FALSE]
    for (v in intersect(active_covs, cat_v)) d[[v]] <- factor(d[[v]])
    covs <- active_covs[vapply(active_covs, function(v) length(unique(d[[v]])) > 1, logical(1))]
    if (!length(covs)) return(NULL)
    dist_m <- dist(mm)
    f <- as.formula(paste("dist_m ~", paste(sprintf("`%s`", covs), collapse = " + ")))
    r <- tryCatch(as.data.frame(vegan::adonis2(f, data = d, by = "margin", permutations = 999)), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    data.table(feature_block = b, covariate = gsub("`", "", rownames(r)), df = r$Df, R2 = r$R2, F = r$F, p_value = r$`Pr(>F)`)[!covariate %in% c("Residual", "Total")]
  }))
  if (nrow(perm)) {
    perm[, note := fifelse(df == 0, "aliased with other covariates (not testable)", "")]
    perm[df == 0, `:=`(F = NA_real_, R2 = NA_real_, p_value = NA_real_)]
    wt(perm, cfg, "table07b_permanova")
  }
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
