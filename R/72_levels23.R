## Binding models, Levels 2 and 3 (step 19), from the Level 1 tables of every modelled sample.
## Level 2: for each read version x model x repeat class x motif family, the per-sample estimates (with their standard
##   errors) are combined by random-effects meta-analysis (metafor::rma, REML):
##     overall   intercept-only: the effect pooled over samples, with tau2 (spread between samples) and I2
##     terms     meta-regression on the model terms (binding_models.terms; default every term step 17 marks estimable):
##               each term's change in the motif effect relative to the reference level, and the omnibus test QM
##   P-values are BH-adjusted within read version x model (overall) and within read version x model x term (terms).
## Level 3: ashr shrinkage of the overall effects (per read version x model) and of each term's effects (per read
##   version x model x term); ranked by the shrunk effect. No significance threshold is applied in the tables.
## Outputs: tables/models_*.tsv and figures/mod01-mod05 (see report/figures.yaml).
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "utils.R")); source(file.path(here, "bin_io.R")); source(file.path(here, "models_lib.R"))
cfg <- args_config(); s <- bm_settings(cfg)
if (!s$enabled) { msg("binding models disabled"); quit(status = 0) }
need_pkgs(c("metafor", "ashr"))
suppressPackageStartupMessages(library(ggplot2))
L1 <- file.path(cfg$work, "models", "level1")
dm <- design_matrix(cfg); terms <- model_terms(cfg, s)
msg("model terms: ", if (nrow(terms)) paste(terms$term, collapse = ", ") else "none (overall effects only)")

rd_all <- function(pat) rbindlist(lapply(Sys.glob(file.path(L1, pat)), function(f) if (file.size(f) > 1) fread(f) else NULL), fill = TRUE)
co <- rd_all("*.coef.tsv"); info <- rd_all("*.info.tsv"); va <- rd_all("*.valid.tsv")
if (!nrow(co)) stop("no Level 1 results in ", L1, ": run step 18 (and its array job) first")
miss <- setdiff(dm$sample_id, unique(co$sample_id))
if (length(miss)) msg("WARNING: no Level 1 results for ", paste(miss, collapse = ", "), "; they are left out of Levels 2-3")
co <- co[sample_id %in% dm$sample_id]
wt(co, cfg, "models_level1_coefficients"); wt(info, cfg, "models_level1_info")

## ---------------------------------------------------------------- Level 2
D <- co[status == "estimated" & is.finite(estimate) & is.finite(se) & se > 0]
D <- merge(D, dm, by = "sample_id")
fit_group <- function(g) {
  k <- nrow(g); base <- list(k = k)
  ov <- tryCatch(metafor::rma(yi = g$estimate, sei = g$se, method = if (k >= 3) "REML" else "FE"), error = function(e) NULL)
  o <- if (is.null(ov)) data.table(k = k, estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
                                   tau2 = NA_real_, I2 = NA_real_, method = "failed") else
    data.table(k = k, estimate = as.numeric(ov$b), se = ov$se, ci_low = ov$ci.lb, ci_high = ov$ci.ub, p = ov$pval,
               tau2 = ov$tau2, I2 = ov$I2, method = ov$method)
  tt <- data.table()
  if (nrow(terms)) {
    use <- terms$term[vapply(terms$term, function(t) length(unique(g[[t]])) > 1, logical(1))]
    dropped <- setdiff(terms$term, use)
    if (length(use) && k > length(use) + 1) {
      X <- as.matrix(g[, ..use])
      mr <- tryCatch(metafor::rma(yi = g$estimate, sei = g$se, mods = X, method = "REML"), error = function(e) NULL)
      if (!is.null(mr)) {
        cf <- data.table(term = sub("^X", "", rownames(mr$beta))[-1], estimate = as.numeric(mr$beta)[-1], se = mr$se[-1],
                         p = mr$pval[-1], qm_p = mr$QMp, tau2_residual = mr$tau2, k = k, status = "estimated")
        cf[, term := use[match(term, make.names(use))] %||% term]
        cf[is.na(term), term := use[seq_len(.N)]]
        tt <- cf
      }
    }
    lack <- setdiff(terms$term, tt$term)
    if (length(lack)) tt <- rbind(tt, data.table(term = lack, estimate = NA_real_, se = NA_real_, p = NA_real_, qm_p = NA_real_,
                                                 tau2_residual = NA_real_, k = k,
                                                 status = ifelse(lack %in% dropped, "no variation among the samples with an estimate",
                                                                 "too few samples for the meta-regression")), fill = TRUE)
  }
  list(o = o, t = tt)
}
keys <- c("arm", "model", "class", "family_id")
G <- split(D, by = keys, keep.by = TRUE)
res <- lapply(G, fit_group)
OV <- rbindlist(lapply(seq_along(G), function(i) cbind(G[[i]][1, ..keys], res[[i]]$o)))
TE <- rbindlist(lapply(seq_along(G), function(i) if (nrow(res[[i]]$t)) cbind(G[[i]][1, ..keys], res[[i]]$t) else NULL), fill = TRUE)
OV[, q := bh(p), by = .(arm, model)]
if (nrow(TE)) TE[, q := bh(p), by = .(arm, model, term)]

## ---------------------------------------------------------------- Level 3 (ashr)
shrink <- function(x) {
  ok <- is.finite(x$estimate) & is.finite(x$se) & x$se > 0
  x[, `:=`(shrunk = NA_real_, shrunk_sd = NA_real_, lfsr = NA_real_)]
  if (sum(ok) >= 3) {
    a <- ashr::ash(x$estimate[ok], x$se[ok], mixcompdist = "normal")
    x[ok, `:=`(shrunk = ashr::get_pm(a), shrunk_sd = ashr::get_psd(a), lfsr = ashr::get_lfsr(a))]
  }
  x
}
OV <- rbindlist(lapply(split(OV, by = c("arm", "model")), shrink))
OV[, rank := frank(-shrunk, na.last = "keep", ties.method = "min"), by = .(arm, model)]
if (nrow(TE)) TE <- rbindlist(lapply(split(TE, by = c("arm", "model", "term")), shrink))
setorder(OV, arm, model, rank, na.last = TRUE)
num4 <- function(x) { for (c in names(x)) if (is.double(x[[c]])) set(x, j = c, value = signif(x[[c]], 4)); x }
wt(num4(copy(OV)), cfg, "models_motif_effects")
if (nrow(TE)) wt(num4(copy(TE)), cfg, "models_condition_effects")

## before and after: pooled effect ignoring repeat class vs within repeat class
BA <- dcast(OV[model %in% c("pooled_unadjusted", "pooled_adjusted")], arm + family_id ~ model,
            value.var = c("estimate", "se", "shrunk", "shrunk_sd", "lfsr"))
if (nrow(BA)) {
  BA[, `:=`(difference = estimate_pooled_unadjusted - estimate_pooled_adjusted,
            share_explained_by_repeat_context = fifelse(abs(estimate_pooled_unadjusted) > 2 * se_pooled_unadjusted,   # only for a clear unadjusted effect
                                                        1 - estimate_pooled_adjusted / estimate_pooled_unadjusted, NA_real_))]
  wt(num4(copy(BA)), cfg, "models_before_after")
}
## read versions: all reads vs unique reads, per class and family
if (all(c("all", "uniq") %in% OV$arm)) {
  AG <- dcast(OV[model == "per_class"], class + family_id ~ arm, value.var = c("estimate", "se", "shrunk"))
  AG[, same_sign := sign(estimate_all) == sign(estimate_uniq)]
  wt(num4(copy(AG)), cfg, "models_read_version_agreement")
}
## held-out arrays, summarised over samples
if (nrow(va)) {
  va <- va[sample_id %in% dm$sample_id]
  wt(num4(copy(va)), cfg, "models_validation_by_sample")
  VS <- va[, .(n_samples = .N, mean_delta_mse = mean(delta_mse, na.rm = TRUE), median_delta_mse = median(delta_mse, na.rm = TRUE),
               samples_improved = sum(delta_mse > 0, na.rm = TRUE), mean_mse_full = mean(mse_full, na.rm = TRUE)),
           by = .(arm, class, family_id)]
  VS[, mean_relative_delta := mean_delta_mse / mean_mse_full]
  wt(num4(copy(VS)), cfg, "models_validation")
}

## ---------------------------------------------------------------- figures
fam_lab <- function(x) x
A1 <- if ("all" %in% OV$arm) "all" else OV$arm[1]
pc <- OV[arm == A1 & model == "per_class" & is.finite(shrunk)]
if (nrow(pc)) {
  topf <- pc[, .(m = max(abs(shrunk))), by = family_id][order(-m)][seq_len(min(30, .N)), family_id]
  pd <- pc[family_id %in% topf][, family_id := factor(family_id, rev(topf))]
  ncl <- uniqueN(pd$class)
  if (ncl <= 6) {
    ## few classes: one panel per class, shrunk effect with its 95% interval
    p1 <- ggplot(pd, aes(shrunk, family_id)) + geom_vline(xintercept = 0, colour = "grey60") +
      geom_errorbar(aes(xmin = shrunk - 1.96 * shrunk_sd, xmax = shrunk + 1.96 * shrunk_sd), width = 0, orientation = "y", colour = "grey40") +
      geom_point(size = 1.6, colour = "#2a78d6") + facet_wrap(~class, nrow = 1) +
      labs(x = "Effect on binding per unit motif term (shrunk, 95% interval)", y = NULL,
           title = sprintf("Motif effects within repeat classes (%s reads)", A1)) + theme_x(9)
    w1 <- max(7, 2.5 + 1.8 * ncl); h1 <- max(4, 1.2 + 0.18 * length(topf))
  } else {
    ## many classes (CenSat + RepeatMasker): one row per family, one column per class; colour = shrunk effect,
    ## dot = lfsr < 0.05 (sign of the effect reliable). One panel per class would exceed ggsave's 50-inch limit.
    lim <- max(abs(pd$shrunk), 1e-6)
    p1 <- ggplot(pd, aes(class, family_id, fill = shrunk)) + geom_tile(colour = "white") +
      geom_point(data = pd[lfsr < 0.05], size = 0.9, colour = "black", show.legend = FALSE) +
      scale_fill_gradient2(low = "#2a78d6", mid = "#f7f7f7", high = "#e34948", limits = c(-lim, lim),
                           name = "Effect per unit\nmotif term\n(shrunk)") +
      labs(x = NULL, y = NULL, caption = "Dot: local false sign rate < 0.05",
           title = sprintf("Motif effects within repeat classes (%s reads)", A1)) + theme_x(9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), legend.position = "right")
    w1 <- max(7, 3 + 0.22 * ncl + text_in(as.character(topf), 9)); h1 <- max(4, 1.8 + 0.18 * length(topf) + text_in(unique(pd$class), 9) * 0.7)
  }
  save_plot(p1, cfg, "mod01_motif_effects", w1, h1)
  fig_note(cfg, "mod01_motif_effects", sprintf("%d families with the largest shrunk effect in any class; %d samples.", length(topf), uniqueN(co$sample_id)))
}
if (exists("BA") && nrow(BA[arm == A1])) {
  bd <- BA[arm == A1 & is.finite(shrunk_pooled_unadjusted)]
  bd <- bd[order(-abs(shrunk_pooled_unadjusted))][seq_len(min(30, .N))]
  bd[, family_id := factor(family_id, rev(family_id))]
  p2 <- ggplot(bd, aes(y = family_id)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_segment(aes(x = shrunk_pooled_unadjusted, xend = shrunk_pooled_adjusted, yend = family_id), colour = "grey50") +
    geom_point(aes(x = shrunk_pooled_unadjusted, colour = "Repeat class ignored"), size = 2) +
    geom_point(aes(x = shrunk_pooled_adjusted, colour = "Within repeat class"), size = 2) +
    scale_colour_manual(values = c("Repeat class ignored" = "#eb6834", "Within repeat class" = "#2a78d6"), name = NULL) +
    labs(x = "Effect on binding per unit motif term (shrunk)", y = NULL, title = "Motif effects before and after the repeat-class comparison") + theme_x(9)
  save_plot(p2, cfg, "mod02_before_after", 7, max(4, 1.4 + 0.18 * nrow(bd)))
  fig_note(cfg, "mod02_before_after", sprintf("%d families; pooled over repeat classes, %s reads.", nrow(bd), A1))
}
if (nrow(TE) && nrow(te <- TE[arm == A1 & model == "pooled_adjusted" & is.finite(shrunk)])) {
  topf <- te[, .(m = max(abs(shrunk))), by = family_id][order(-m)][seq_len(min(30, .N)), family_id]
  te <- te[family_id %in% topf][, family_id := factor(family_id, rev(topf))]
  lim <- max(abs(te$shrunk), 1e-6)
  p3 <- ggplot(te, aes(term, family_id, fill = shrunk)) + geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f7f7f7", high = "#e34948", limits = c(-lim, lim), name = "Change in\nmotif effect") +
    labs(x = NULL, y = NULL, title = "Condition effects on motif effects (pooled over repeat classes, within class)") + theme_x(9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), legend.position = "right")
  save_plot(p3, cfg, "mod03_condition_effects", max(6, 3 + 0.6 * uniqueN(te$term)), max(4, 1.6 + 0.18 * length(topf)))
  fig_note(cfg, "mod03_condition_effects", sprintf("%d terms x %d families (largest shrunk change).", uniqueN(te$term), length(topf)))
}
if (exists("VS") && nrow(VS)) {
  vd <- VS[arm == A1]
  p4 <- ggplot(vd, aes(class, family_id, fill = mean_relative_delta)) + geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f7f7f7", high = "#1baf7a", name = "Held-out MSE\nreduction\n(relative)") +
    labs(x = NULL, y = NULL, title = "Held-out arrays: prediction improvement from each motif family") + theme_x(9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), legend.position = "right")
  save_plot(p4, cfg, "mod04_validation", max(5, 3 + 0.22 * uniqueN(vd$class) + text_in(unique(vd$family_id), 9)),
            max(4, 1.4 + 0.16 * uniqueN(vd$family_id) + 0.7 * text_in(unique(vd$class), 9)))
  fig_note(cfg, "mod04_validation", sprintf("%d CenSat classes with at least %d arrays; mean over samples.", uniqueN(vd$class), s$validation_folds))
}
if (exists("AG") && nrow(ag5 <- AG[is.finite(estimate_all) & is.finite(estimate_uniq)])) {
  p5 <- ggplot(ag5, aes(estimate_all, estimate_uniq)) +
    geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    geom_abline(linetype = 2, colour = "grey50") + geom_point(alpha = 0.6, colour = "#2a78d6") + facet_wrap(~class, ncol = nc5 <- min(uniqueN(ag5$class), 6)) +
    labs(x = "Effect, all reads", y = "Effect, uniquely placed reads", title = "Motif effects: all reads vs unique reads") + theme_x(9)
  save_plot(p5, cfg, "mod05_read_versions", max(6, 1 + 2 * nc5), max(4, 1 + 2 * ceiling(uniqueN(ag5$class) / nc5)))
  fig_note(cfg, "mod05_read_versions", sprintf("%d family x class effects with both read versions.", nrow(AG[is.finite(estimate_all) & is.finite(estimate_uniq)])))
}
msg(sprintf("Levels 2-3: %d family x class x model x read-version groups; %d condition-effect estimates", nrow(OV), if (nrow(TE)) TE[status == "estimated", .N] else 0L))
