## Binding models: similarity between samples (step 20).
## Binding (log2 IP/input) of every modelled sample in one read version (binding_models.similarity_arm) on the bins that
## pass the filters in every sample, separately for CenSat bins and all other bins (at most similarity_max_bins of each,
## sampled with binding_models.seed). For each pair: Pearson correlation r_ij. Each sample's ceiling is its split-half
## reliability: the correlation of binding from the two random halves of its reads, stepped up to full depth
## (Spearman-Brown: 2r / (1 + r)). Scaled similarity = r_ij / sqrt(ceiling_i x ceiling_j).
## Similarity by condition, on the distance 1 - scaled similarity:
##   PERMANOVA (vegan::adonis2, binding_models.permutations, marginal test of each model term)
##   one-way ICC per covariate: (MSB - MSW) / (MSB + (n0 - 1) MSW) from the same sums of squares, and its R2
## Outputs: tables/models_similarity_pairs.tsv, models_similarity_ceilings.tsv, models_similarity_permanova.tsv,
##          models_similarity_icc.tsv; figures/mod06_similarity_<group>.png, mod07_similarity_by_condition.png
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "utils.R")); source(file.path(here, "bin_io.R")); source(file.path(here, "models_lib.R"))
cfg <- args_config(); s <- bm_settings(cfg)
if (!s$enabled) { msg("binding models disabled"); quit(status = 0) }
need_pkgs("vegan"); suppressPackageStartupMessages({ library(ggplot2); library(pheatmap) })
B <- bin_dir(cfg); P <- model_libraries(cfg, B); dm <- design_matrix(cfg); terms <- model_terms(cfg, s)
arm <- s$similarity_arm; if (!arm %in% c("all", "uniq")) stop("binding_models.similarity_arm must be all or uniq")
A <- load_model_annot(cfg, B, s)
keep_seq <- A$keep_seq; is_censat <- A$is_censat; rm(A); invisible(gc())
## Memory: genome-wide vectors are never held for all samples at once (11 samples x 3 vectors x tens of millions of
## bins exhausts a node). Pass 1 builds the bins passing in every sample; the bins are then chosen (subsample with the
## seed); pass 2 recomputes each sample and keeps only the chosen bins. The signal is read twice.
sample_vals <- function(i) {
  b <- sample_binding(P[i], arm, s)
  if (is.null(b)) return(NULL)
  ip <- load_signal(P$ip_dir[i], cols = c("h1", "h2")); ct <- load_signal(P$ctl_dir[i], cols = c("h1", "h2"))
  Li <- attr(ip, "library"); Lc <- attr(ct, "library")
  y1 <- binding(ip$h1, ct$h1, Li$h1, Lc$h1); y2 <- binding(ip$h2, ct$h2, Li$h2, Lc$h2)
  ## split-half ceiling: each input half must pass the input filter on its own. Selecting bins on the full input
  ## (half 1 + half 2) makes the halves anti-correlated when coverage is low (selection on a sum), biasing the ceiling down.
  hk <- function(x, n) { v <- x / (n / 1e6); t <- as.numeric(quantile(v[v > 0], s$input_min_quantile, na.rm = TRUE)); !is.na(v) & v > 0 & v >= t }
  hok <- hk(ct$h1, Lc$h1) & hk(ct$h2, Lc$h2)
  y1[!hok] <- NA; y2[!hok] <- NA
  list(y = b$y, y1 = y1, y2 = y2, ok = b$input_ok & is.finite(b$y))
}
ok <- keep_seq; sm <- character()
for (i in seq_len(nrow(P))) {
  v <- sample_vals(i)
  if (is.null(v)) { msg("WARNING: ", P$sample_id[i], ": no reads in read version ", arm, "; left out"); next }
  ok <- ok & v$ok; sm <- c(sm, P$sample_id[i]); rm(v); invisible(gc())
}
if (length(sm) < 3) stop("similarity needs at least 3 samples with binding in read version ", arm)
msg(sprintf("similarity: %d samples, %s reads; %d bins pass the filters in every sample", length(sm), arm, sum(ok)))
set.seed(s$seed)
IDX <- list()
for (g in c("censat", "other")) {
  idx <- which(ok & (is_censat == (g == "censat")))
  if (length(idx) < 100) { msg("similarity, ", g, " bins: only ", length(idx), " bins; skipped"); next }
  if (length(idx) > s$similarity_max_bins) idx <- sort(sample(idx, s$similarity_max_bins))
  IDX[[g]] <- idx
}
rm(ok, keep_seq, is_censat); invisible(gc())
Y <- setNames(lapply(names(IDX), function(g) list()), names(IDX)); H1 <- Y; H2 <- Y
for (i in which(P$sample_id %in% sm)) {
  v <- sample_vals(i)
  for (g in names(IDX)) { k <- IDX[[g]]; Y[[g]][[P$sample_id[i]]] <- v$y[k]; H1[[g]][[P$sample_id[i]]] <- v$y1[k]; H2[[g]][[P$sample_id[i]]] <- v$y2[k] }
  rm(v); invisible(gc())
}
pairs <- list(); ceil <- list(); perm <- list(); icc <- list(); mats <- list(); used <- list()
meta <- samples_dt(cfg)[match(sm, sample_id)]
for (g in names(IDX)) {
  idx <- IDX[[g]]
  M <- sapply(sm, function(x) Y[[g]][[x]])
  R <- cor(M)
  h <- vapply(sm, function(x) { a <- H1[[g]][[x]]; b2 <- H2[[g]][[x]]; k <- is.finite(a) & is.finite(b2); if (sum(k) < 100) NA_real_ else cor(a[k], b2[k]) }, numeric(1))
  cl <- 2 * h / (1 + h)
  ## a sample whose two read halves do not correlate positively has no usable ceiling: it is kept in the r table but
  ## left out of scaled similarity and the tests on it (listed here and in models_similarity_ceilings$status)
  rel <- is.finite(cl) & cl > 0
  ceil[[g]] <- data.table(group = g, sample_id = sm, split_half_r = h, ceiling = cl, n_bins = length(idx),
                          status = fifelse(rel, "used", "left out of scaled similarity: split-half r <= 0 or fewer than 100 bins with both input halves passing"))
  if (any(!rel)) msg("WARNING: similarity, ", g, " bins: split-half r <= 0 (no reliable signal) for ", paste(sm[!rel], collapse = ", "),
                     "; left out of scaled similarity, PERMANOVA and ICC")
  S <- R / sqrt(outer(pmax(cl, 0), pmax(cl, 0))); S[!rel, ] <- NA; S[, !rel] <- NA; diag(S)[rel] <- 1
  ut <- which(upper.tri(R), arr.ind = TRUE)
  pairs[[g]] <- data.table(group = g, sample_a = sm[ut[, 1]], sample_b = sm[ut[, 2]], r = R[ut], scaled = S[ut], n_bins = length(idx))
  sg <- sm[rel]
  if (length(sg) < 3) { msg("similarity, ", g, " bins: fewer than 3 samples with a usable ceiling; no scaled similarity tests"); next }
  S <- S[rel, rel, drop = FALSE]; mats[[g]] <- S; used[[g]] <- sg
  Dm <- as.dist(1 - pmin(S, 1))
  tcols <- terms$term; X <- as.data.frame(dm[match(sg, sample_id), ..tcols])
  use <- names(X)[vapply(X, function(v) length(unique(v)) > 1, logical(1))]
  if (length(use) && length(sg) > length(use) + 1) {
    names(X) <- make.names(names(X))
    f <- stats::as.formula(paste("Dm ~", paste(make.names(use), collapse = " + ")))
    ad <- as.data.frame(vegan::adonis2(f, data = X, permutations = s$permutations, by = "margin"))
    ad$term <- rownames(ad)
    perm[[g]] <- data.table(group = g, term = ad$term, df = ad$Df, sum_sq = ad$SumOfSqs, R2 = ad$R2, F = ad$F, p = ad$`Pr(>F)`,
                            permutations = s$permutations, n_samples = length(sg))[!term %in% c("Residual", "Total")][, term := use[match(term, make.names(use))]]
  }
  mg <- meta[match(sg, sample_id)]
  for (v in unique(terms$covariate)) {
    if (!v %in% names(meta)) next
    grp <- factor(mg[[v]]); if (nlevels(grp) < 2) next
    ad1 <- as.data.frame(vegan::adonis2(Dm ~ grp, permutations = s$permutations))
    ssb <- ad1$SumOfSqs[1]; ssw <- ad1$SumOfSqs[2]; N <- length(sg); g_ <- nlevels(grp)
    msb <- ssb / (g_ - 1); msw <- ssw / (N - g_); n0 <- (N - sum(table(grp)^2) / N) / (g_ - 1)
    icc[[paste(g, v)]] <- data.table(group = g, covariate = v, levels = g_, R2 = ad1$R2[1], icc = (msb - msw) / (msb + (n0 - 1) * msw),
                                     p = ad1$`Pr(>F)`[1], permutations = s$permutations, n_samples = N)
  }
}
r4 <- function(x) { x <- rbindlist(x, fill = TRUE); for (c in names(x)) if (is.double(x[[c]])) set(x, j = c, value = signif(x[[c]], 4)); x }
wt(r4(pairs), cfg, "models_similarity_pairs"); wt(r4(ceil), cfg, "models_similarity_ceilings")
if (length(perm)) wt(r4(perm), cfg, "models_similarity_permanova")
if (length(icc)) wt(r4(icc), cfg, "models_similarity_icc")
for (g in names(mats)) {
  S <- mats[[g]]; sg <- used[[g]]; dimnames(S) <- list(short_id(sg), short_id(sg))
  ann <- sample_annotation(cfg, meta[match(sg, sample_id)], sg); rownames(ann) <- short_id(sg)
  w <- max(6, 3 + 0.3 * length(sg) + text_in(rownames(S), 8)); h <- max(5, 2 + 0.3 * length(sg) + text_in(rownames(S), 8))
  save_pheatmap(function() pheatmap(S, annotation_col = ann, fontsize = 8, color = colorRampPalette(c("#f7f7f7", "#2a78d6", "#0c447c"))(50),
                                    main = sprintf("Scaled similarity of binding, %s bins", if (g == "censat") "CenSat" else "other"),
                                    cellwidth = 0.3 * 72, cellheight = 0.3 * 72),
                cfg, paste0("mod06_similarity_", g), width = w, height = h)
  fig_note(cfg, paste0("mod06_similarity_", g), sprintf("%d samples, %d bins, %s reads.%s", length(sg), ceil[[g]]$n_bins[1], arm, if (length(sg) < length(sm)) sprintf(" %d sample(s) with split-half r <= 0 left out.", length(sm) - length(sg)) else ""))
}
if (length(icc)) {
  ic <- rbindlist(icc)
  p7 <- ggplot(ic, aes(covariate, R2, fill = group)) + geom_col(position = position_dodge(0.8), width = 0.7) +
    scale_fill_manual(values = c(censat = "#eb6834", other = "#2a78d6"), labels = c(censat = "CenSat bins", other = "Other bins"), name = NULL) +
    labs(x = NULL, y = "Share of between-sample variation (one-way R2)", title = "Similarity of binding explained by each condition") + theme_x(9)
  save_plot(p7, cfg, "mod07_similarity_by_condition", max(5, 2 + 0.8 * uniqueN(ic$covariate)), 4.5)
  fig_note(cfg, "mod07_similarity_by_condition", sprintf("%d covariates; %d permutations.", uniqueN(ic$covariate), s$permutations))
}
