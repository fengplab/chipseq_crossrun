## Regression models: which motifs are driven by repeats, and what varies with sample metadata.
## Input: work/models/ (step 16), work/motifs/families.tsv, work/scan/*.family_sea.tsv, work/reads/reads_in_labels.tsv
##
## A. Motif family selection    table08_family_selection
## B. Confounding (with run)    table10_covariate_associations, table10b_design_aliasing, table10c_variance_partition,
##                              fig17_design_confounding
## C. Peak-level mixed logistic model, one per family
##                              table09_motif_classification, table09b_peak_model_coefficients,
##                              table09c_peak_model_covariates, table09d_peak_model_excluded_classes, fig18_motif_classification
## D. Sample-level quasi-binomial / quasi-Poisson models against metadata
##                              table11_covariate_tests, table11b_covariate_coefficients, table11c_stratified_tests,
##                              fig19_covariate_effects_<feature type>
## E. Time-point contrasts within each series (bootstrap over peaks)
##                              table12_timepoint_contrasts, fig20_timepoint_contrasts_<series>
## Run ID is used only in the confounding analysis, the restricted permutations and the run-sensitivity refits.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
suppressPackageStartupMessages({ library(lme4); library(patchwork) })
cfg <- args_config()
mc <- cfg$models %||% list()
opt <- function(k, d) { v <- mc[[k]]; if (is.null(v)) d else v }
if (!isTRUE(opt("enabled", TRUE))) { msg("models disabled (models.enabled: false)"); quit(status = 0) }
set.seed(as.integer(opt("seed", 1)))
alpha <- as.numeric(opt("alpha", 0.05))
max_rows <- as.integer(opt("heatmap_max_rows", 60))   # rows per heatmap part (Figures 18-19)
threads <- max(1L, as.integer(cfg$threads %||% 1))
MD <- file.path(cfg$work, "models")
clear_figures(cfg, "^fig(17|18|19|20)_")
for (t in list.files(cfg$tables, "^table(08|09|10|11|12)", full.names = TRUE)) unlink(t)
summary_note <- function(fig, text) fig_note(cfg, paste0(fig, "__summary"), text)
for (f in c("fig17", "fig18", "fig19", "fig20")) fig_note(cfg, paste0(f, "__summary"), "")
drop_empty_notes <- function() {
  f <- file.path(cfg$figures, "figure_notes.tsv")
  if (file.exists(f)) { n <- fread(f, sep = "\t", quote = ""); fwrite(n[!is.na(note) & note != ""], f, sep = "\t", quote = FALSE) }
}

S <- samples_dt(cfg)
S <- S[peaks_file != ""]
## run / batch variable: optional (batch_var: null = not used). It is never part of the confounding analysis; when set, it
## is used only to restrict permutations and for the run-sensitivity refits.
bv <- cfg$batch_var
if (is.null(bv) || !nzchar(bv) || !bv %in% names(S)) bv <- NULL
sv <- cfg$series_var; tv <- cfg$time_var
cat_v <- intersect(unlist(cfg$categorical_covariates), names(S))
cov_lab <- function(v) {
  if (v == "time_bin") return(paste0(cov_lab(tv), " (binned)"))
  if (!is.null(bv) && v == bv) return("Pipeline run")
  l <- cfg$covariate_labels[[v]]
  if (!is.null(l)) return(l)
  gsub("_", " ", v)
}
rdgz <- function(f) fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t", quote = "", na.strings = c("", "NA"))

## ======================= 0. load windows =========================================================
wf <- list.files(file.path(MD, "win"), "\\.windows\\.tsv\\.gz$", full.names = TRUE)
if (!length(wf)) stop("no window tables in work/models/win; run step 16 first")
W <- rbindlist(lapply(wf, rdgz), fill = TRUE)
W <- W[sample_id %in% S$sample_id]
S <- S[sample_id %in% W$sample_id]
W[, peak := as.integer(set == "peak")]
LIB <- rd(file.path(MD, "libsize.tsv"))
W <- merge(W, LIB, by = "sample_id", all.x = TRUE)
cmin <- as.numeric(opt("class_min_fraction", 0.5))
W[, rclass := fifelse(censat_class != "none" & censat_frac >= cmin, paste0("censat:", censat_class),
               fifelse(rm_class != "none" & rm_frac >= cmin, paste0("repeatmasker:", rm_class), "none"))]
class_label <- function(x) {
  fs <- sub(":.*$", "", x); lb <- sub("^[^:]*:", "", x)
  ifelse(x == "none", "Outside repeats", ifelse(fs == "censat", pretty_censat(lb), paste0("RepeatMasker ", lb)))
}
class_name <- function(x) { x <- as.character(x); out <- rep(NA_character_, length(x)); ok <- !is.na(x); out[ok] <- class_label(x[ok]); out }
## IP/input ratio of a window (log2, library-size normalised) and its excess over the sample's baseline:
## the baseline is a low quantile of the ratio among the sample's own peaks, i.e. roughly the ratio needed to be called
W[, log2_ratio := log2(((ip_reads + 0.5) / ip_total_reads) / ((ctl_reads + 0.5) / ctl_total_reads))]
bq <- as.numeric(opt("ratio_baseline_quantile", 0.10))
W[, ratio_baseline := if (any(peak == 1 & !is.na(log2_ratio))) quantile(log2_ratio[peak == 1], bq, na.rm = TRUE) else NA_real_, by = sample_id]
W[, ratio_excess := fifelse(peak == 1, log2_ratio - ratio_baseline, 0)]
msg(sprintf("windows: %d peaks (%d in the peak-level model) and %d background windows from %d samples",
            W[peak == 1, .N], W[peak == 1 & in_model == 1, .N], W[peak == 0, .N], uniqueN(W$sample_id)))

## ======================= A. motif family selection ===============================================
sf <- list.files(file.path(MD, "win"), "\\.sites\\.tsv\\.gz$", full.names = TRUE)
wkey <- W[, .(sample_id, window_id, peak, in_model, width, rclass)]
setkey(wkey, sample_id, window_id)
## pass 1: per sample x family aggregates for every family
AG <- rbindlist(lapply(sf, function(f) {
  d <- rdgz(f); if (!nrow(d)) return(NULL)
  d <- merge(d, wkey, by = c("sample_id", "window_id"))
  d[, .(pk_present = sum(peak == 1 & in_model == 1), bg_present = sum(peak == 0),
        sites_pk = sum(n_sites[peak == 1]), tandem_pk = sum(n_tandem_sites[peak == 1]), sites_bg = sum(n_sites[peak == 0])),
    by = .(sample_id, family_id)]
}))
TOT <- W[, .(n_pk_model = sum(peak == 1 & in_model == 1), n_bg = sum(peak == 0), n_pk = sum(peak == 1),
             kb_pk = sum(width[peak == 1]) / 1000, kb_bg = sum(width[peak == 0]) / 1000), by = sample_id]
fam <- rd(file.path(cfg$work, "motifs", "families.tsv"))
feat <- rd(file.path(cfg$work, "motifs", "family_rep_features.tsv"))
FS <- fam[, .(n_members = .N, n_samples_discovered = uniqueN(sample_id),
              min_member_evalue = suppressWarnings(min(as.numeric(evalue), na.rm = TRUE)),
              min_member_sea_pvalue = suppressWarnings(min(as.numeric(sea_pvalue), na.rm = TRUE))), by = family_id]
FS[!is.finite(min_member_evalue), min_member_evalue := NA_real_][!is.finite(min_member_sea_pvalue), min_member_sea_pvalue := NA_real_]
FS[, evalue_size_corrected := min_member_evalue * n_members]
## combined SEA p-value across samples (Cauchy combination: valid under arbitrary dependence between samples)
cauchy <- function(p) {
  p <- p[!is.na(p)]; if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-300), 1 - 1e-12)
  t <- mean(ifelse(p < 1e-15, 1 / (p * pi), tan((0.5 - p) * pi)))
  if (t > 1e15) 1 / (t * pi) else 0.5 - atan(t) / pi
}
fsea <- rbindlist(lapply(list.files(file.path(cfg$work, "scan"), "\\.family_sea\\.tsv$", full.names = TRUE), function(f) {
  l <- readLines(f); l <- l[!grepl("^#", l) & nzchar(l)]
  if (length(l) < 2) return(NULL)
  d <- fread(text = paste(l, collapse = "\n"))
  setnames(d, tolower(gsub("%", "_pct", gsub("[^A-Za-z0-9%]+", "_", names(d)))))
  p <- as.numeric(d$pvalue)
  if ("log_pvalue" %in% names(d)) p <- ifelse(!is.na(p) & p == 0 & as.numeric(d$log_pvalue) < 0, exp(as.numeric(d$log_pvalue)), p)
  data.table(sample_id = sub("\\.family_sea\\.tsv$", "", basename(f)), family_id = d$id, p = p)
}), fill = TRUE)
if (nrow(fsea)) FS <- merge(FS, fsea[sample_id %in% S$sample_id, .(sea_combined_pvalue = cauchy(p), n_samples_sea = sum(!is.na(p))), by = family_id],
                            by = "family_id", all.x = TRUE) else FS[, `:=`(sea_combined_pvalue = NA_real_, n_samples_sea = 0L)]
## pooled peaks-vs-background test, stratified by sample (Cochran-Mantel-Haenszel)
cmh <- function(a, b, n1, n0) {     # a/n1 = present/total in peaks, b/n0 in background, one element per sample
  a <- as.numeric(a); b <- as.numeric(b); n1 <- as.numeric(n1); n0 <- as.numeric(n0)
  ok <- n1 > 0 & n0 > 0 & (a + b) > 0 & (a + b) < (n1 + n0)
  a <- a[ok]; b <- b[ok]; n1 <- n1[ok]; n0 <- n0[ok]
  if (!length(a)) return(list(or = NA_real_, p = NA_real_, k = 0L))
  N <- n1 + n0
  or <- sum(a * (n0 - b) / N) / sum(b * (n1 - a) / N)
  if (length(a) == 1) return(list(or = or, p = fisher.test(matrix(c(a, n1 - a, b, n0 - b), 2))$p.value, k = 1L))
  e <- n1 * (a + b) / N; v <- n1 * n0 * (a + b) * (N - a - b) / (N^2 * (N - 1))
  chi <- (abs(sum(a - e)) - 0.5)^2 / sum(v)
  list(or = or, p = pchisq(chi, 1, lower.tail = FALSE), k = length(a))
}
grid <- CJ(sample_id = TOT$sample_id, family_id = FS$family_id)
AGf <- merge(grid, AG, by = c("sample_id", "family_id"), all.x = TRUE)
for (v in c("pk_present", "bg_present", "sites_pk", "tandem_pk", "sites_bg")) set(AGf, which(is.na(AGf[[v]])), v, 0)
AGf <- merge(AGf, TOT, by = "sample_id")
PO <- AGf[, { r <- cmh(pk_present, bg_present, n_pk_model, n_bg)
              .(pooled_odds_ratio = r$or, pooled_pvalue = r$p, pooled_n_samples = r$k,
                pct_peaks_with_motif = 100 * sum(pk_present) / sum(n_pk_model), pct_background_with_motif = 100 * sum(bg_present) / sum(n_bg)) },
          by = family_id]
PO[, pooled_qvalue := bh(pooled_pvalue)]
FS <- merge(FS, PO, by = "family_id", all.x = TRUE)
FS <- merge(feat[, .(family_id, pwm_consensus)], FS, by = "family_id")
stat <- opt("family_stat", "sea_pvalue"); cutoff <- mc$family_cutoff
FS[, selection_stat := stat]
FS[, selection_value := switch(stat, sea_pvalue = sea_combined_pvalue, evalue = evalue_size_corrected,
                               recurrence = as.numeric(n_samples_discovered), pooled_q = pooled_qvalue,
                               stop("models.family_stat must be sea_pvalue, evalue, recurrence or pooled_q"))]
FS[, passes_cutoff := !is.na(selection_value) & (if (is.null(cutoff)) TRUE else if (stat == "recurrence") selection_value >= cutoff else selection_value <= cutoff)]
setorderv(FS, "selection_value", order = if (stat == "recurrence") -1L else 1L, na.last = TRUE)
maxf <- if ("max_families" %in% names(mc) && is.null(mc$max_families)) Inf else as.numeric(opt("max_families", 50))
FS[, selected := passes_cutoff & cumsum(passes_cutoff) <= maxf]
wt(FS, cfg, "table08_family_selection")
sel <- FS[selected == TRUE, family_id]
msg(sprintf("family selection by %s%s: %d of %d families pass, %d enter the models (max_families = %s)", stat,
            if (is.null(cutoff)) " (no cutoff)" else sprintf(" %s %g", if (stat == "recurrence") ">=" else "<=", cutoff),
            sum(FS$passes_cutoff), nrow(FS), length(sel), if (is.finite(maxf)) format(maxf) else "no limit"))
fam_label <- function(f) paste(f, FS$pwm_consensus[match(f, FS$family_id)])
## pass 2: sites of the selected families
SI <- rbindlist(lapply(sf, function(f) { d <- rdgz(f); d[family_id %in% sel] }))
## sanity check: motif sites must join to the windows. If they do not (e.g. FIMO reported chromosome names instead of
## window names), every family would look absent and the models would be meaningless.
if (length(sel)) {
  joined <- merge(SI[, .(sample_id, window_id, family_id)], wkey[, .(sample_id, window_id)], by = c("sample_id", "window_id"))
  fam_ok <- uniqueN(joined$family_id)
  if (!nrow(joined) || fam_ok < 0.5 * length(sel))
    stop(sprintf(paste("only %d of %d selected motif families have any site in any window (%d site rows, %d joined to a window).",
                       "The motif sites do not match the window IDs; rerun steps 03-16 with the current bundle",
                       "(FIMO >= 5.5.4 must not parse FASTA headers as genomic coordinates)."),
                 fam_ok, length(sel), nrow(SI), nrow(joined)))
}

## ======================= B. design and confounding ==============================================
## ordered time categories shared across series
tb <- sort(as.numeric(unlist(opt("time_bins", c(4, 12)))))
if (!is.null(tv) && tv %in% names(S) && any(!is.na(S[[tv]]))) {
  labs <- c(paste0("<=", tb[1]), if (length(tb) > 1) paste0(">", head(tb, -1), "-", tb[-1]), paste0(">", tail(tb, 1)))
  S[, time_bin := as.character(cut(get(tv), c(-Inf, tb, Inf), labels = labs))]
  S[is.na(time_bin), time_bin := "no time"]
  S[, time_bin := factor(time_bin, levels = c(labs, "no time")[c(labs, "no time") %in% time_bin])]
}
covs <- c(cat_v, if ("time_bin" %in% names(S)) "time_bin")
for (v in cat_v) S[, (v) := factor(fifelse(is.na(get(v)) | get(v) == "", "NA", as.character(get(v))))]
covs <- covs[vapply(covs, function(v) nlevels(droplevels(S[[v]])) > 1, logical(1))]
if (!is.null(bv)) S[, (bv) := factor(get(bv))]

## design matrix of the covariates (treatment contrasts); aliased columns are dropped and reported
mm_of <- function(d, v) { m <- model.matrix(~ x, data.frame(x = droplevels(d[[v]])))[, -1, drop = FALSE]
                          colnames(m) <- paste0(v, "__", levels(droplevels(d[[v]]))[-1]); m }
build_design <- function(d, vars) {
  cols <- list(); owner <- character()
  for (v in vars) { if (nlevels(droplevels(d[[v]])) < 2) next
                    m <- mm_of(d, v); cols[[v]] <- m; owner <- c(owner, rep(v, ncol(m))) }
  if (!length(cols)) return(list(Z = matrix(0, nrow(d), 0), owner = character(), dropped = NULL))
  Z <- do.call(cbind, cols); keep <- rep(TRUE, ncol(Z)); dropped <- list()
  for (j in seq_len(ncol(Z))) {                       # keep a column only if it adds rank to the intercept + earlier columns
    A <- cbind(1, Z[, which(keep & seq_len(ncol(Z)) < j), drop = FALSE])
    if (qr(cbind(A, Z[, j]))$rank == qr(A)$rank) {
      keep[j] <- FALSE
      co <- qr.coef(qr(A), Z[, j]); co <- co[-1]; co <- co[!is.na(co) & abs(co) > 1e-8]
      dropped[[length(dropped) + 1]] <- data.table(term = colnames(Z)[j], covariate = owner[j],
                                                    aliased_with = if (length(co)) paste(names(co), collapse = " + ") else "intercept (constant)")
    }
  }
  list(Z = Z[, keep, drop = FALSE], owner = owner[keep], dropped = rbindlist(dropped))
}
DES <- build_design(S, covs); Z <- DES$Z; rownames(Z) <- S$sample_id
term_label <- function(term) {
  v <- sub("__.*$", "", term); l <- sub("^.*?__", "", term)
  sprintf("%s: %s vs %s", vapply(v, cov_lab, ""), l, vapply(v, function(x) levels(droplevels(S[[x]]))[1], ""))
}
## pairwise association (Cramer's V), including run
cramers <- function(a, b) { tb <- table(droplevels(a), droplevels(b)); k <- min(dim(tb)) - 1
  if (k < 1) return(NA_real_); sqrt(unname(suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)) / (sum(tb) * k)) }
av <- covs
assoc <- rbindlist(if (length(av) >= 2) combn(av, 2, simplify = FALSE, FUN = function(p) {
  v <- cramers(S[[p[1]]], S[[p[2]]])
  data.table(variable_a = p[1], variable_b = p[2], label_a = cov_lab(p[1]), label_b = cov_lab(p[2]), cramers_v = v, n = nrow(S),
             note = fifelse(!is.na(v) & v >= 0.99, "aliased: effects cannot be separated",
                    fifelse(!is.na(v) & v >= 0.7, "strongly associated: interpret effects with caution", "")))
}) else list())
if (nrow(assoc)) wt(assoc[order(-cramers_v)], cfg, "table10_covariate_associations")
## aliasing table: one row per covariate level
varies_within <- function(v, block) !is.null(block) && block %in% names(S) && any(S[, uniqueN(get(v)), by = c(block)]$V1 > 1)
fb <- opt("block_fallback", "cell_type"); if (!fb %in% names(S)) fb <- NULL
al <- rbindlist(lapply(covs, function(v) {
  lv <- levels(droplevels(S[[v]]))
  d <- data.table(covariate = v, covariate_label = cov_lab(v), level = lv, n_samples = as.integer(table(droplevels(S[[v]]))[lv]),
                  term = c(NA, paste0(v, "__", lv[-1])))
  d[, status := fifelse(is.na(term), "reference level", fifelse(term %in% colnames(Z), "estimable", "dropped (aliased)"))]
  d[, aliased_with := if (!is.null(DES$dropped) && nrow(DES$dropped)) DES$dropped$aliased_with[match(term, DES$dropped$term)] else NA_character_]
  d[, permutation_block := if (varies_within(v, bv)) bv else if (!is.null(fb) && fb != v && varies_within(v, fb)) fb else "none (unrestricted)"]
  d[, term := NULL][]
}))
if (nrow(al)) wt(al, cfg, "table10b_design_aliasing")
if (!is.null(DES$dropped) && nrow(DES$dropped))
  msg("aliased covariate levels dropped from the models: ", paste(sprintf("%s (= %s)", DES$dropped$term, DES$dropped$aliased_with), collapse = "; "))
msg(sprintf("design: %d samples, %d covariates (%s), %d estimable terms", nrow(S), length(covs), paste(covs, collapse = ", "), ncol(Z)))

## ======================= C. peak-level mixed logistic model ====================================
mcell <- as.integer(opt("min_cell_count", 5)); confined <- as.numeric(opt("confined_fraction", 0.95))
use_region <- isTRUE(opt("region_weighting", TRUE))
MW <- W[in_model == 1 & (peak == 1 | bg_draw == "1" | bg_draw == 1) & !is.na(log2_ratio) & !is.na(ratio_baseline) & !is.na(gc)]
no_ratio <- setdiff(S$sample_id, unique(MW$sample_id))
if (length(no_ratio)) msg("samples without IP/control read counts are left out of the peak-level model: ", paste(no_ratio, collapse = ", "))
## The same region can be a peak in several samples, and its sequence (hence motif presence) is then the same each time.
## To stop such regions being counted as independent evidence, every peak is weighted by 1 / (number of modelled
## samples in which its merged region is a peak), so each region contributes a total weight of 1. (A random effect per
## region is not used: motif presence is nearly fixed within a region, which makes the peak coefficient describe a
## "typical region" and can reverse its sign relative to the overall enrichment.)
MW[, region := fifelse(peak == 0 | is.na(region_id), window_id, as.character(region_id))]
MW[, w_region := if (use_region) 1 / .N else 1, by = .(region, peak)]
MW[, `:=`(log2_width = log2(width) - mean(log2(width)), gc10 = (gc - mean(gc)) * 10)]
MW[, sample_id := factor(sample_id)]
## covariate x peak columns; the covariates are centred over samples, so the 'peak' term is the peak effect of the average sample
Zc <- scale(Z, center = TRUE, scale = FALSE)
Zp <- Zc[as.character(MW$sample_id), , drop = FALSE] * MW$peak
if (ncol(Zp)) colnames(Zp) <- paste0("cov_", make.names(colnames(Z)))
cov_terms <- data.table(col = colnames(Zp), term = colnames(Z), covariate = DES$owner)
if (!is.null(bv)) MW[, run := S[[bv]][match(as.character(sample_id), S$sample_id)]]
wald <- function(b, V, idx) {      # joint Wald chi-square test of a set of coefficients
  idx <- idx[!is.na(b[idx])]; if (!length(idx)) return(c(chisq = NA_real_, df = 0, p = NA_real_))
  x <- tryCatch(as.numeric(t(b[idx]) %*% solve(V[idx, idx, drop = FALSE], b[idx])), error = function(e) NA_real_)
  c(chisq = x, df = length(idx), p = pchisq(x, length(idx), lower.tail = FALSE))
}
fit_family <- function(fid) {
  d <- copy(MW)
  hit <- SI[family_id == fid, .(sample_id, window_id, n_sites)]
  d[, y := as.integer(paste(sample_id, window_id) %in% paste(hit$sample_id, hit$window_id))]
  base <- list(family_id = fid, n_windows = nrow(d), n_peaks = sum(d$peak), n_motif_peaks = sum(d$y[d$peak == 1]),
               n_motif_background = sum(d$y[d$peak == 0]),
               n_regions = uniqueN(d$region[d$peak == 1]), n_motif_regions = uniqueN(d$region[d$peak == 1 & d$y == 1]))
  ## which repeat classes can be modelled: every cell of (peak/background) x (motif present/absent) needs >= mcell windows
  cells <- d[, .(pk_present = sum(peak == 1 & y == 1), pk_absent = sum(peak == 1 & y == 0),
                 bg_present = sum(peak == 0 & y == 1), bg_absent = sum(peak == 0 & y == 0)), by = rclass]
  cells[, ok := pmin(pk_present, pk_absent, bg_present, bg_absent) >= mcell]
  in_rep <- d[y == 1, mean(rclass != "none")]
  base$pct_motif_windows_in_repeats <- 100 * in_rep
  excl <- cells[rclass != "none" & !ok, .(family_id = fid, repeat_class = rclass, pk_present, pk_absent, bg_present, bg_absent,
                                           reason = fifelse(pk_present + bg_present < mcell, "motif (almost) absent in this class",
                                                    fifelse(pk_absent + bg_absent < mcell, "motif in (almost) every window of this class",
                                                    fifelse(pk_present + pk_absent < 2 * mcell, "too few peaks in this class", "too few background windows in this class"))))]
  none_ok <- cells[rclass == "none", isTRUE(ok)]
  if (!isTRUE(none_ok)) {
    top <- cells[rclass != "none"][order(-(pk_present + bg_present))][1]
    cls <- if (sum(d$y) >= mcell && !is.na(in_rep) && in_rep >= confined) "confined to repeats" else "not modelled"
    return(list(summary = c(base, list(classification = cls, top_repeat_class = if (cls == "confined to repeats") top$rclass else NA_character_,
                model_note = if (cls == "confined to repeats")
                  sprintf("%.0f%% of motif-containing windows lie in repeats, too few elsewhere to fit a model", 100 * in_rep)
                else "too few windows with (or without) the motif outside repeats to fit a model")), excl = excl))
  }
  keep_cl <- cells[rclass != "none" & ok, rclass]
  d <- d[rclass %in% c("none", keep_cl)]
  zp <- Zp[MW$rclass %in% c("none", keep_cl), , drop = FALSE]
  zp <- zp[, colSums(abs(zp)) > 0 & apply(zp, 2, function(x) length(unique(x[d$peak == 1])) > 1), drop = FALSE]
  ccol <- character(); icol <- character()
  for (k in seq_along(keep_cl)) {
    set(d, j = paste0("cls_", k), value = as.integer(d$rclass == keep_cl[k]))
    set(d, j = paste0("int_", k), value = as.integer(d$rclass == keep_cl[k]) * d$peak)
    ccol <- c(ccol, paste0("cls_", k)); icol <- c(icol, paste0("int_", k))
  }
  if (ncol(zp)) d <- cbind(d, as.data.table(zp))
  fixed <- c("peak", "ratio_excess", ccol, icol, "log2_width", "gc10", colnames(zp))
  ## random effects: sample intercept and sample-specific peak effect, so that covariate x peak terms are judged against
  ## between-sample variation in peak enrichment
  re_opts <- c("(1 + peak | sample_id)", "(1 | sample_id)")
  fit <- NULL; note <- ""
  for (re in re_opts) {
    fit <- tryCatch(suppressWarnings(suppressMessages(glmer(as.formula(paste("y ~", paste(fixed, collapse = " + "), "+", re)), data = d, family = binomial, weights = w_region, nAGQ = 0,
                           control = glmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")))), error = function(e) NULL)
    if (!is.null(fit)) { if (re != re_opts[1]) note <- paste0("random effects simplified to ", re); break }
  }
  if (is.null(fit)) return(list(summary = c(base, list(classification = "not modelled", model_note = "model could not be fitted")), excl = excl))
  b <- fixef(fit); V <- as.matrix(vcov(fit)); se <- sqrt(diag(V))
  co <- data.table(family_id = fid, term = names(b), estimate = unname(b), std_error = unname(se))
  co[, p_value := 2 * pnorm(-abs(estimate / std_error))]
  co[, term_label := fifelse(term == "(Intercept)", "Intercept (background window outside repeats)",
                     fifelse(term == "peak", "Peak vs background: outside repeats, at the baseline IP/input ratio, average sample",
                     fifelse(term == "ratio_excess", "Per doubling of IP/input above the sample's baseline (peaks)",
                     fifelse(term == "log2_width", "Per doubling of window width",
                     fifelse(term == "gc10", "Per 10 percentage points of GC", term)))))]
  for (k in seq_along(keep_cl)) {
    co[term == paste0("cls_", k), term_label := paste0(class_label(keep_cl[k]), ": motif frequency vs outside repeats (sequence effect)")]
    co[term == paste0("int_", k), term_label := paste0(class_label(keep_cl[k]), ": extra peak enrichment inside this class")]
  }
  if (nrow(cov_terms)) for (i in seq_len(nrow(cov_terms)))
    co[term == cov_terms$col[i], term_label := paste0("Peak enrichment, ", term_label(cov_terms$term[i]))]
  co[, term_group := fifelse(grepl("^cls_", term), "repeat class", fifelse(grepl("^int_", term), "repeat class x peak",
                     fifelse(grepl("^cov_", term), "covariate x peak", fifelse(term == "peak", "peak", "adjustment"))))]
  co[, repeat_class := NA_character_]
  for (k in seq_along(keep_cl)) co[term %in% paste0(c("cls_", "int_"), k), repeat_class := keep_cl[k]]
  wr <- wald(b, V, ccol); wi <- wald(b, V, icol)
  cv <- rbindlist(lapply(unique(cov_terms$covariate), function(v) {
    idx <- intersect(cov_terms[covariate == v, col], names(b)); if (!length(idx)) return(NULL)
    w <- wald(b, V, idx)
    data.table(family_id = fid, covariate = v, covariate_label = cov_lab(v), chisq = w[["chisq"]], df = w[["df"]], p_value = w[["p"]],
               largest_effect = cov_terms$term[match(idx[which.max(abs(b[idx]))], cov_terms$col)], largest_log_odds = b[idx][which.max(abs(b[idx]))])
  }))
  ## run sensitivity: let the peak effect also vary by run
  if (!is.null(cv) && nrow(cv) && !is.null(bv) && isTRUE(opt("run_sensitivity", TRUE)) && uniqueN(d$run) > 1) {
    f2 <- tryCatch(suppressWarnings(suppressMessages(glmer(as.formula(paste("y ~", paste(fixed, collapse = " + "), "+", re_opts[1], "+ (0 + peak | run)")), data = d, family = binomial, weights = w_region, nAGQ = 0,
                          control = glmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")))), error = function(e) NULL)
    if (!is.null(f2)) { b2 <- fixef(f2); V2 <- as.matrix(vcov(f2))
      cv[, p_value_run_adjusted := vapply(covariate, function(v) wald(b2, V2, intersect(cov_terms[covariate == v, col], names(b2)))[["p"]], 0)] }
  }
  top_c <- if (length(ccol)) keep_cl[which.max(abs(b[ccol]))] else NA_character_
  top_i <- if (length(icol)) keep_cl[which.max(abs(b[icol]))] else NA_character_
  list(summary = c(base, list(n_windows_modelled = nrow(d), n_repeat_classes = length(keep_cl),
                              peak_log_odds = b[["peak"]], peak_std_error = se[["peak"]], peak_p = co[term == "peak", p_value],
                              repeat_chisq = wr[["chisq"]], repeat_df = wr[["df"]], repeat_p = wr[["p"]],
                              interaction_chisq = wi[["chisq"]], interaction_df = wi[["df"]], interaction_p = wi[["p"]],
                              ratio_slope = b[["ratio_excess"]], ratio_p = co[term == "ratio_excess", p_value],
                              top_repeat_class = top_c, top_repeat_log_odds = if (length(ccol)) b[ccol][which.max(abs(b[ccol]))] else NA_real_,
                              top_interaction_class = top_i, top_interaction_log_odds = if (length(icol)) b[icol][which.max(abs(b[icol]))] else NA_real_,
                              sample_sd_intercept = tryCatch(attr(VarCorr(fit)$sample_id, "stddev")[[1]], error = function(e) NA_real_),
                              model_note = note)),
       coef = co, cov = cv, excl = excl)
}
if (length(sel) && nrow(MW)) {
  msg(sprintf("peak-level mixed logistic models: %d families x %d windows (%d samples)", length(sel), nrow(MW), uniqueN(MW$sample_id)))
  res <- parallel::mclapply(sel, function(f) tryCatch(fit_family(f), error = function(e)
    list(summary = list(family_id = f, classification = "not modelled", model_note = paste("error:", conditionMessage(e))))), mc.cores = min(threads, length(sel)))
  T9 <- rbindlist(lapply(res, function(r) as.data.table(r$summary)), fill = TRUE)
  for (v in c("peak_log_odds", "peak_std_error", "peak_p", "repeat_chisq", "repeat_df", "repeat_p", "interaction_chisq", "interaction_df",
              "interaction_p", "ratio_slope", "ratio_p", "top_repeat_class", "top_repeat_log_odds", "top_interaction_class",
              "top_interaction_log_odds", "n_windows_modelled", "n_repeat_classes", "sample_sd_intercept", "pct_motif_windows_in_repeats"))
    if (!v %in% names(T9)) T9[, (v) := if (grepl("class$", v)) NA_character_ else NA_real_]
  T9[, `:=`(peak_q = bh(peak_p), repeat_q = bh(repeat_p), interaction_q = bh(interaction_p), ratio_q = bh(ratio_p))]
  sig <- function(q) !is.na(q) & q < alpha
  T9[is.na(classification), classification := fifelse(sig(interaction_q), "repeat-specific binding",
                                               fifelse(sig(peak_q) & peak_log_odds > 0, "universal",
                                               fifelse(sig(peak_q) & peak_log_odds < 0, "depleted in peaks",
                                               fifelse(sig(repeat_q), "repeat-driven", "no association"))))]
  T9[classification == "universal" & sig(repeat_q), model_note := trimws(paste(model_note, "motif frequency also differs between repeat classes"))]
  T9[classification == "no association" & n_repeat_classes == 0, model_note := trimws(paste(model_note, "no repeat class had enough windows to test"))]
  T9[, top_repeat_class := as.character(top_repeat_class)][, top_interaction_class := as.character(top_interaction_class)]
  T9[, top_repeat_class_name := class_name(top_repeat_class)]
  T9[, top_interaction_class_name := class_name(top_interaction_class)]
  rl <- rd(file.path(cfg$tables, "table02_repeat_likeness.tsv"))
  T9[, pwm_consensus := FS$pwm_consensus[match(family_id, FS$family_id)]]
  T9[, repeat_likeness_flag := if (!is.null(rl) && "flag" %in% names(rl)) rl$flag[match(family_id, rl$family_id)] else NA_character_]
  T9 <- merge(T9, FS[, .(family_id, pooled_odds_ratio, pooled_qvalue)], by = "family_id", all.x = TRUE)
  setcolorder(T9, c("family_id", "pwm_consensus", "classification", "repeat_likeness_flag"))
  wt(T9[order(classification, pmin(peak_q, repeat_q, interaction_q, na.rm = TRUE))], cfg, "table09_motif_classification")
  CO <- rbindlist(lapply(res, `[[`, "coef"), fill = TRUE)
  if (nrow(CO)) { CO[, repeat_class := as.character(repeat_class)][, repeat_class_name := class_name(repeat_class)]
    CO[, q_value := bh(p_value), by = term_group]
    wt(CO[, .(family_id, pwm_consensus = FS$pwm_consensus[match(family_id, FS$family_id)], term_group, term_label, repeat_class, repeat_class_name,
              estimate, std_error, odds_ratio = exp(estimate), p_value, q_value)], cfg, "table09b_peak_model_coefficients") }
  CV <- rbindlist(lapply(res, `[[`, "cov"), fill = TRUE)
  if (nrow(CV)) { CV[, q_value := bh(p_value)]
    if ("p_value_run_adjusted" %in% names(CV)) CV[, q_value_run_adjusted := bh(p_value_run_adjusted)]
    CV[, largest_effect := term_label(largest_effect)]
    CV[, pwm_consensus := FS$pwm_consensus[match(family_id, FS$family_id)]]
    setcolorder(CV, c("family_id", "pwm_consensus"))
    wt(CV[order(p_value)], cfg, "table09c_peak_model_covariates") }
  EX <- rbindlist(lapply(res, `[[`, "excl"), fill = TRUE)
  if (nrow(EX)) { EX[, repeat_class_name := class_label(repeat_class)]
    wt(EX[, .(family_id, pwm_consensus = FS$pwm_consensus[match(family_id, FS$family_id)], repeat_class, repeat_class_name,
              pk_present, pk_absent, bg_present, bg_absent, reason)], cfg, "table09d_peak_model_excluded_classes") }
  msg("motif classification: ", paste(sprintf("%s = %d", names(table(T9$classification)), as.integer(table(T9$classification))), collapse = ", "))

  ## ---- Fig 18: what explains motif presence. Rows: for every column, all families significant there, or that column's
  ## 5 smallest p-values if fewer are significant (union over columns); split into parts of at most heatmap_max_rows rows.
  clear_figures(cfg, "^fig18_")
  pd <- if (nrow(CO)) CO[term_group %in% c("peak", "repeat class", "repeat class x peak")] else CO
  if (nrow(pd)) {
    pd[, col := fifelse(term_group == "peak", "In peaks\n(outside repeats)", stringr_wrap(class_label(repeat_class), 22))]
    pd[, col_key := paste(term_group, col)]
    ## a tile is starred when its own q-value is below alpha and, for repeat-class terms, the family's joint test of that
    ## group of terms (the test used for the classification) is also significant, so stars and labels agree
    jq <- T9[, .(family_id, repeat_q, interaction_q)]
    pd <- merge(pd, jq, by = "family_id", all.x = TRUE)
    pd[, is_sig := !is.na(q_value) & q_value < alpha &
         (term_group == "peak" | (term_group == "repeat class" & !is.na(repeat_q) & repeat_q < alpha) |
          (term_group == "repeat class x peak" & !is.na(interaction_q) & interaction_q < alpha))]
    showf <- select_heatmap_rows(pd, "family_id", "col_key", "p_value", "is_sig", 5)
    pd <- pd[family_id %in% showf]
    pd[, col := factor(col, levels = c("In peaks\n(outside repeats)", setdiff(sort(unique(col)), "In peaks\n(outside repeats)")))]
    pd[, panel := factor(fifelse(term_group == "repeat class", "Sequence effect:\nmotif frequency in the class", "Binding effect:\npeak vs background"),
                         levels = c("Sequence effect:\nmotif frequency in the class", "Binding effect:\npeak vs background"))]
    cls_of <- function(f) T9$classification[match(f, T9$family_id)]
    ## order: by classification, then strongest evidence first
    showf <- showf[order(cls_of(showf), match(showf, showf))]
    pd[, fam := paste0(fam_label(family_id), "  [", cls_of(family_id), "]")]
    lim <- max(2, min(8, max(abs(pd$estimate), na.rm = TRUE)))
    pages <- paginate(showf, max_rows)
    for (i in seq_along(pages)) {
      pf <- pages[[i]]; dp <- pd[family_id %in% pf]
      lv <- paste0(fam_label(pf), "  [", cls_of(pf), "]")
      dp[, fam := factor(fam, levels = rev(lv))]
      p18 <- ggplot(dp, aes(col, fam)) + geom_tile(aes(fill = pmax(pmin(estimate, lim), -lim)), colour = "grey85") +
        geom_point(data = dp[is_sig == TRUE], shape = 8, size = 1.6) +
        scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#d6402a", midpoint = 0, limits = c(-lim, lim), name = "log odds ratio") +
        scale_y_discrete(drop = FALSE) + facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
        labs(x = NULL, y = NULL, title = "Peak-level mixed logistic models: what explains motif presence",
             subtitle = if (length(pages) > 1) sprintf("Part %d of %d", i, length(pages)) else NULL) +
        theme_x(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(), legend.key.width = grid::unit(1.2, "cm"))
      nm <- page_name("fig18_motif_classification", i, length(pages))
      save_plot(p18, cfg, nm, max(8, 4.5 + 0.55 * uniqueN(pd$col)), 2.8 + 0.3 * length(pf))
      fig_note(cfg, nm, sprintf(
        "%d families in total (part %d of %d shows %d): for each column, every family with q < %g there, or the 5 smallest p-values when fewer are significant. Classification over all %d selected families: %s. Families confined to repeats or not modelled have no coefficients and are listed in table09. Star: q < %g for the term and, for repeat-class terms, for the family's joint test of those terms (as used for the classification).",
        length(showf), i, length(pages), length(pf), alpha, nrow(T9),
        paste(sprintf("%s %d", names(table(T9$classification)), as.integer(table(T9$classification))), collapse = ", "), alpha))
    }
  } else summary_note("fig18", "No motif family could be modelled at the peak level (see table09_motif_classification for the reasons).")
} else {
  summary_note("fig18", if (!length(sel)) sprintf("No motif family passed the selection cutoff (models.family_stat = %s); see table08_family_selection.", stat)
                        else "No sample has IP and control read counts, which the peak-level model needs.")
}

## ======================= D. sample-level models against metadata ================================
## Every response is a count out of a total (quasi-binomial) or a count per exposure (quasi-Poisson), so samples with
## more peaks/reads/sites weigh more, and the dispersion is estimated from the data.
R <- list()
add_resp <- function(type, id, label, d, family) if (nrow(d)) R[[length(R) + 1]] <<- list(type = type, id = id, label = label, d = d, family = family)
safe_logit <- function(x, n) qlogis((x + 0.5) / (n + 1))
## 1. IP reads per repeat class (control fraction as offset)
rr <- rd(file.path(cfg$work, "reads", "reads_in_labels.tsv"))
if (!is.null(rr) && nrow(rr)) {
  ip <- merge(S[bam_file != "", .(sample_id, bam = bam_file, control_bam)], rr, by = "bam")
  ip <- merge(ip, rr[, .(control_bam = bam, feature_set, label, c_reads = reads, c_total = total_reads)], by = c("control_bam", "feature_set", "label"))
  for (k in unique(paste(ip$feature_set, ip$label, sep = ":"))) {
    d <- ip[paste(feature_set, label, sep = ":") == k]
    if (sum(d$reads) < 50) next
    add_resp("reads_class", k, paste0("IP reads in ", pretty_label(d$feature_set[1], d$label[1]), " (vs control)"),
             d[, .(sample_id, y = reads, n = total_reads, off = safe_logit(c_reads, c_total))], "quasibinomial")
  }
}
## 2. IP reads in peaks (control reads in the same peaks as offset)
pk <- W[peak == 1 & !is.na(ip_reads) & !is.na(ctl_reads), .(y = sum(ip_reads), n = ip_total_reads[1], cy = sum(ctl_reads), cn = ctl_total_reads[1]), by = sample_id]
add_resp("reads_peaks", "peaks", "IP reads in peaks (vs control reads in the same peaks)",
         pk[, .(sample_id, y = pmin(y, n), n, off = safe_logit(pmin(cy, cn), cn))], "quasibinomial")
## 3. peak share per repeat class
for (k in setdiff(unique(W[peak == 1, rclass]), "none")) {
  d <- W[peak == 1, .(y = sum(rclass == k), n = .N), by = sample_id]
  if (sum(d$y) < 10) next
  add_resp("peak_share", k, paste0("Share of peaks in ", class_label(k)), d[, .(sample_id, y, n, off = 0)], "quasibinomial")
}
## 4-6. motif responses for the selected families
if (length(sel)) {
  for (f in sel) {
    a <- AGf[family_id == f]
    add_resp("tandem", f, paste0("Tandem array score of ", fam_label(f)), a[sites_pk > 0, .(sample_id, y = tandem_pk, n = sites_pk, off = 0)], "quasibinomial")
    add_resp("density", f, paste0("Motif density per kb in peaks of ", fam_label(f), " (vs background)"),
             a[kb_pk > 0 & kb_bg > 0, .(sample_id, y = sites_pk, n = NA_real_, off = log(kb_pk) + log((sites_bg + 0.5) / kb_bg))], "quasipoisson")
  }
  PS <- merge(SI, W[peak == 1, .(sample_id, window_id, rclass)], by = c("sample_id", "window_id"))
  cs <- PS[rclass != "none", .(y = sum(n_sites)), by = .(sample_id, family_id, rclass)]
  ck <- W[peak == 1 & rclass != "none", .(kb = sum(width) / 1000), by = .(sample_id, rclass)]
  for (k in unique(ck$rclass)) for (f in sel) {
    d <- merge(ck[rclass == k], cs[rclass == k & family_id == f, .(sample_id, y)], by = "sample_id", all.x = TRUE)
    d[is.na(y), y := 0]
    if (sum(d$y) < 10) next
    add_resp("density_class", paste(f, k, sep = "|"), paste0("Density of ", fam_label(f), " per kb of peaks in ", class_label(k)),
             d[, .(sample_id, y, n = NA_real_, off = log(kb))], "quasipoisson")
  }
}
TYPE_LAB <- c(reads_class = "IP vs control reads per repeat class", reads_peaks = "IP vs control reads in peaks",
              peak_share = "Share of peaks per repeat class", tandem = "Motif tandem array score",
              density = "Motif density per kb in peaks", density_class = "Motif density per kb inside each repeat class")
msg(sprintf("sample-level responses: %s", paste(sprintf("%s %d", names(table(vapply(R, `[[`, "", "type"))), as.integer(table(vapply(R, `[[`, "", "type")))), collapse = ", ")))

min_n <- as.integer(opt("min_samples", 6)); nperm <- as.integer(opt("n_permutations", 499)); perm_p <- as.numeric(opt("perm_max_p", 0.05))
glm_fit <- function(d, X, family) {
  X <- cbind(`(Intercept)` = 1, X)
  suppressWarnings(tryCatch(
    if (family == "quasibinomial") glm.fit(X, cbind(d$y, d$n - d$y), family = quasibinomial(), offset = d$off)
    else glm.fit(X, d$y, family = quasipoisson(), offset = d$off), error = function(e) NULL))
}
f_test <- function(d, X, cols, family) {       # F test of the columns `cols` of X given the others
  f1 <- glm_fit(d, X, family); f0 <- glm_fit(d, X[, setdiff(colnames(X), cols), drop = FALSE], family)
  if (is.null(f1) || is.null(f0) || !isTRUE(f1$converged) || f1$df.residual < 1) return(c(F = NA_real_, df = NA_real_, p = NA_real_, disp = NA_real_))
  disp <- sum(f1$weights * f1$residuals^2) / f1$df.residual
  dfn <- f0$df.residual - f1$df.residual
  if (dfn < 1 || !is.finite(disp) || disp <= 0) return(c(F = NA_real_, df = dfn, p = NA_real_, disp = disp))
  Fv <- max(0, (f0$deviance - f1$deviance) / dfn) / disp
  c(F = Fv, df = dfn, p = pf(Fv, dfn, f1$df.residual, lower.tail = FALSE), disp = disp)
}
usable <- function(X, d) {                     # drop design columns that are constant or aliased in this subset of samples
  if (!ncol(X)) return(X)
  keep <- logical(ncol(X))
  for (j in seq_len(ncol(X))) { A <- cbind(1, X[, which(keep), drop = FALSE]); keep[j] <- qr(cbind(A, X[, j]))$rank > qr(A)$rank }
  X[, keep, drop = FALSE]
}
perm_block <- function(v, ids) {               # blocks for restricted permutation of covariate v
  s <- S[match(ids, sample_id)]
  for (b in c(bv, if (!is.null(fb) && fb != v) fb)) if (!is.null(b) && any(s[, uniqueN(get(v)), by = c(b)]$V1 > 1)) return(list(name = b, g = as.character(s[[b]])))
  list(name = "none (unrestricted)", g = rep("all", length(ids)))
}
fit_response <- function(r) {
  d <- r$d[sample_id %in% rownames(Z)]
  d <- d[is.finite(off) & !is.na(y)]
  if (nrow(d) < min_n || uniqueN(d$y) < 2) return(NULL)
  X <- usable(Z[d$sample_id, , drop = FALSE], d)
  if (!ncol(X) || nrow(d) - ncol(X) - 1 < 2) return(NULL)
  own <- DES$owner[match(colnames(X), colnames(Z))]
  full <- glm_fit(d, X, r$family)
  if (is.null(full) || !isTRUE(full$converged)) return(NULL)
  disp <- sum(full$weights * full$residuals^2) / full$df.residual
  se <- sqrt(diag(chol2inv(full$qr$qr[seq_len(full$rank), seq_len(full$rank), drop = FALSE])) * disp)
  co <- data.table(term = colnames(X), covariate = own, estimate = full$coefficients[-1], std_error = se[-1])
  co[, p_value := 2 * pt(-abs(estimate / std_error), full$df.residual)]
  Rm <- if (!is.null(bv)) mm_of(S[match(d$sample_id, sample_id)], bv) else NULL
  tests <- rbindlist(lapply(unique(own), function(v) {
    cols <- colnames(X)[own == v]; ft <- f_test(d, X, cols, r$family)
    out <- data.table(covariate = v, F = ft[["F"]], df = ft[["df"]], df_residual = full$df.residual, p_value = ft[["p"]], dispersion = disp,
                      p_value_within_run = NA_real_, run_note = "", perm_p_value = NA_real_, perm_block = NA_character_)
    ## run sensitivity: the same test with run as a blocking factor (only possible if the covariate varies within runs)
    if (!is.null(Rm) && ncol(Rm)) {
      others <- usable(cbind(Rm, X[, setdiff(colnames(X), cols), drop = FALSE]), d)
      Xr <- cbind(others, X[, cols, drop = FALSE])
      gain <- qr(cbind(1, Xr))$rank - qr(cbind(1, others))$rank
      if (gain < length(cols)) out[, run_note := if (gain == 0) "constant within runs: cannot be separated from run" else "partly aliased with run"]
      if (gain == length(cols) && nrow(d) - ncol(Xr) - 1 >= 1) out[, p_value_within_run := f_test(d, Xr, cols, r$family)[["p"]]]
      else if (gain == length(cols)) out[, run_note := "too few samples to add run"]
    }
    if (!is.na(ft[["p"]]) && ft[["p"]] <= perm_p && nperm > 0) {
      pb <- perm_block(v, d$sample_id); out[, perm_block := pb$name]
      Fp <- vapply(seq_len(nperm), function(i) {
        idx <- seq_len(nrow(d)); for (g in unique(pb$g)) { w <- which(pb$g == g); if (length(w) > 1) idx[w] <- sample(w) }
        Xp <- X; Xp[, cols] <- X[idx, cols]
        f_test(d, Xp, cols, r$family)[["F"]] }, 0)
      out[, perm_p_value := (1 + sum(Fp >= ft[["F"]] - 1e-12, na.rm = TRUE)) / (1 + sum(!is.na(Fp)))]
    }
    out
  }))
  ## the same covariate alone, within each level of the stratifying variable
  sb <- opt("stratify_by", "cell_type"); strat <- NULL
  if (!is.null(sb) && sb %in% names(S)) strat <- rbindlist(lapply(levels(droplevels(S[[sb]])), function(L) {
    ids <- intersect(d$sample_id, S[get(sb) == L, sample_id]); if (length(ids) < 4) return(NULL)
    ds <- d[sample_id %in% ids]; if (uniqueN(ds$y) < 2) return(NULL)
    rbindlist(lapply(setdiff(covs, sb), function(v) {
      ss <- S[match(ds$sample_id, sample_id)]; if (nlevels(droplevels(ss[[v]])) < 2) return(NULL)
      Xs <- mm_of(ss, v); if (nrow(ds) - ncol(Xs) - 1 < 2) return(NULL)
      ft <- f_test(ds, Xs, colnames(Xs), r$family)
      data.table(stratum = L, covariate = v, n = nrow(ds), levels = paste(levels(droplevels(ss[[v]])), collapse = ", "), F = ft[["F"]], df = ft[["df"]], p_value = ft[["p"]])
    }))
  }))
  hdr <- data.table(feature_type = r$type, response_id = r$id, response = r$label, model = r$family, n = nrow(d))
  list(tests = cbind(hdr, tests), coef = cbind(hdr, co), strat = if (!is.null(strat) && nrow(strat)) cbind(hdr[, !"n"], strat) else NULL)
}
if (length(R) && ncol(Z)) {
  GR <- parallel::mclapply(R, function(r) tryCatch(fit_response(r), error = function(e) NULL), mc.cores = threads)
  T11 <- rbindlist(lapply(GR, `[[`, "tests"), fill = TRUE)
  if (nrow(T11)) {
    T11[, q_value := bh(p_value), by = feature_type]
    T11[, `:=`(feature_type_label = TYPE_LAB[feature_type], covariate_label = vapply(covariate, cov_lab, ""))]
    if (is.null(bv)) T11[, c("p_value_within_run", "run_note") := NULL]
    setcolorder(T11, intersect(c("feature_type", "feature_type_label", "response_id", "response", "covariate", "covariate_label", "model", "n", "F", "df", "df_residual",
                       "p_value", "q_value", "perm_p_value", "perm_block", "p_value_within_run", "run_note", "dispersion"), names(T11)))
    wt(T11[order(feature_type, p_value)], cfg, "table11_covariate_tests")
    C11 <- rbindlist(lapply(GR, `[[`, "coef"), fill = TRUE)
    C11[, `:=`(term_label = term_label(term), ratio = exp(estimate))]
    C11 <- merge(C11, T11[, .(feature_type, response_id, covariate, covariate_p = p_value, covariate_q = q_value)], by = c("feature_type", "response_id", "covariate"), all.x = TRUE)
    wt(C11[order(feature_type, covariate_p), .(feature_type, response_id, response, model, n, covariate, term_label, estimate, std_error, ratio, p_value, covariate_p, covariate_q)],
       cfg, "table11b_covariate_coefficients")
    ST <- rbindlist(lapply(GR, `[[`, "strat"), fill = TRUE)
    if (nrow(ST)) { ST[, q_value := bh(p_value), by = feature_type]; ST[, covariate_label := vapply(covariate, cov_lab, "")]
      wt(ST[order(feature_type, p_value)], cfg, "table11c_stratified_tests") }
    msg(sprintf("sample-level models: %d responses, %d covariate tests, %d with q < %g", uniqueN(T11[, paste(feature_type, response_id)]), nrow(T11), T11[!is.na(q_value) & q_value < alpha, .N], alpha))
    ## ---- Fig 19: covariate effects, one sub-figure per feature type. Rows: for every covariate level, all responses
    ## significant there, or that level's 5 smallest p-values if fewer are significant; parts of at most heatmap_max_rows rows.
    clear_figures(cfg, "^fig19_")
    for (ft in intersect(names(TYPE_LAB), unique(C11$feature_type))) {
      cc <- C11[feature_type == ft & !is.na(estimate)]
      if (!nrow(cc)) next
      cc[, is_sig := !is.na(covariate_q) & covariate_q < alpha & p_value < 0.05]
      top <- select_heatmap_rows(cc, "response_id", "term_label", "p_value", "is_sig", 5)
      cc <- cc[response_id %in% top]
      cc[, tl := factor(stringr_wrap(term_label, 24))]
      lim <- max(1, min(6, max(abs(cc$estimate) / log(2), na.rm = TRUE)))
      pages <- paginate(top, max_rows)
      for (i in seq_along(pages)) {
        pr <- pages[[i]]; dp <- cc[response_id %in% pr]
        dp[, resp := factor(stringr_wrap(response, 48), levels = rev(unique(stringr_wrap(cc$response[match(pr, cc$response_id)], 48))))]
        p19 <- ggplot(dp, aes(tl, resp)) + geom_tile(aes(fill = pmax(pmin(estimate / log(2), lim), -lim)), colour = "grey85") +
          geom_point(data = dp[is_sig == TRUE], shape = 8, size = 1.6) +
          scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#d6402a", midpoint = 0, limits = c(-lim, lim),
                               name = if (dp$model[1] == "quasibinomial") "log2 odds ratio" else "log2 rate ratio") +
          scale_x_discrete(drop = FALSE) +
          labs(x = NULL, y = NULL, title = TYPE_LAB[[ft]], subtitle = if (length(pages) > 1) sprintf("Part %d of %d", i, length(pages)) else NULL) + theme_x(9) +
          theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(), legend.key.width = grid::unit(1.2, "cm"))
        nm <- page_name(paste0("fig19_covariate_effects_", ft), i, length(pages))
        save_plot(p19, cfg, nm, max(8, 4.6 + 0.75 * nlevels(cc$tl)), 3 + 0.42 * length(pr))
        fig_note(cfg, nm, sprintf(
          "%d of %d responses in total (part %d of %d shows %d): for each covariate level, every response marked significant there, or the 5 smallest p-values when fewer are. %d covariate tests in this feature type have q < %g. Star: the covariate's F test has q < %g and this level differs from the reference (p < 0.05).",
          length(top), uniqueN(C11[feature_type == ft, response_id]), i, length(pages), length(pr),
          T11[feature_type == ft & !is.na(q_value) & q_value < alpha, .N], alpha, alpha))
      }
    }
  } else summary_note("fig19", "No sample-level model could be fitted (too few samples per response, or no residual degrees of freedom).")
} else summary_note("fig19", "No sample-level model could be fitted (no estimable covariate terms or no responses).")

## ---- variance partition: each covariate alone vs jointly with the other covariates -------------------------
vp <- NULL
Ym <- NULL
if (length(R)) {
  Y <- rbindlist(lapply(R, function(r) r$d[, .(sample_id, id = paste(r$type, r$id, sep = "::"),
        v = if (r$family == "quasibinomial") qlogis((y + 0.5) / (n + 1)) - off else log(y + 0.5) - off)]))
  Ym <- dcast(Y, sample_id ~ id, value.var = "v"); ids <- Ym$sample_id
  Ym <- as.matrix(Ym[, -1]); rownames(Ym) <- ids
  Ym <- Ym[intersect(S$sample_id, ids), , drop = FALSE]
  Ym <- Ym[, colSums(is.na(Ym)) == 0 & apply(Ym, 2, function(x) isTRUE(sd(x) > 0)), drop = FALSE]
  if (ncol(Ym) >= 2 && nrow(Ym) >= 6) Ym <- scale(Ym) else Ym <- NULL
}
if (!is.null(Ym) && length(covs)) {
  Ss <- S[match(rownames(Ym), sample_id)]
  adjr2 <- function(X) { if (is.null(X) || !ncol(X)) return(0); X <- usable(X, NULL); n <- nrow(Ym); p <- qr(cbind(1, X))$rank - 1
    if (n - p - 1 < 1) return(NA_real_)
    fit <- lm.fit(cbind(1, X), Ym); r2 <- 1 - sum(fit$residuals^2) / sum(scale(Ym, scale = FALSE)^2); 1 - (1 - r2) * (n - 1) / (n - p - 1) }
  vp <- rbindlist(lapply(covs, function(v) {
    if (nlevels(droplevels(Ss[[v]])) < 2) return(NULL)
    X <- mm_of(Ss, v); oth <- setdiff(covs, v); oth <- oth[vapply(oth, function(o) nlevels(droplevels(Ss[[o]])) > 1, logical(1))]
    O <- if (length(oth)) do.call(cbind, lapply(oth, function(o) mm_of(Ss, o))) else NULL
    a_x <- adjr2(X); a_xo <- adjr2(cbind(X, O)); a_o <- adjr2(O)
    data.table(covariate = v, covariate_label = cov_lab(v), n_samples = nrow(Ym), n_features = ncol(Ym),
               total = a_x,                              # everything the covariate could explain on its own
               unique = a_xo - a_o,                      # what it explains beyond all other covariates
               shared_with_other_covariates = a_x - (a_xo - a_o),
               all_covariates = a_xo)
  }))
  vp[, note := fifelse(is.na(unique), "too few samples for the number of terms",
               fifelse(unique < 0.01 & total > 0.05, "explains variance only jointly with other covariates: cannot be attributed", ""))]
  wt(vp, cfg, "table10c_variance_partition")
}
## ---- Fig 17: confounding
if (nrow(assoc)) {
  lv <- vapply(av, cov_lab, "")
  a2 <- rbind(assoc[, .(a = label_a, b = label_b, v = cramers_v)], assoc[, .(a = label_b, b = label_a, v = cramers_v)], data.table(a = lv, b = lv, v = 1))
  a2[, `:=`(a = factor(a, levels = lv), b = factor(b, levels = rev(lv)))]
  pA <- ggplot(a2, aes(a, b, fill = v)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", v)), size = 3) +
    scale_fill_gradient(low = "white", high = "#d6402a", limits = c(0, 1), name = "Cramer's V") +
    labs(x = NULL, y = NULL, title = "How strongly the covariates are tied") + theme_x(9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
  if (!is.null(vp) && nrow(vp[!is.na(unique)])) {
    v2 <- melt(vp[!is.na(total)], id.vars = "covariate_label", measure.vars = c("unique", "shared_with_other_covariates"), variable.name = "part", value.name = "r2")
    v2 <- v2[!is.na(r2)]; v2[, r2 := pmax(r2, 0)]
    v2[, part := factor(part, levels = c("shared_with_other_covariates", "unique"),
                        labels = c("Shared with other covariates", "Unique to this covariate"))]
    v2[, covariate_label := factor(covariate_label, levels = rev(lv))]
    pB <- ggplot(v2, aes(r2, covariate_label, fill = part)) + geom_col(width = 0.7) +
      scale_fill_manual(values = c("Unique to this covariate" = "#2a78d6", "Shared with other covariates" = "#b9b8b2"), name = NULL) +
      labs(x = "Share of variance in the sample-level responses (adjusted R2)", y = NULL, title = "Variance each covariate can explain") +
      theme_x(9) + guides(fill = guide_legend(nrow = 2, reverse = TRUE))
    p17 <- pA + pB + plot_layout(widths = c(1, 1.1))
  } else p17 <- pA
  save_plot(p17, cfg, "fig17_design_confounding", if (inherits(p17, "patchwork")) 13 else 7, 2.8 + 0.55 * length(lv))
  nd <- if (!is.null(DES$dropped) && nrow(DES$dropped)) paste0("Levels dropped because they are aliased: ", paste(term_label(DES$dropped$term), collapse = "; "), ". ") else "No covariate level is fully aliased with the others. "
  nr <- ""
  fig_note(cfg, "fig17_design_confounding", paste0(sprintf("%d samples. ", nrow(S)), nd, nr,
           if (is.null(vp)) " Variance partition not available (too few samples or responses)." else ""))
}

## ======================= E. time-point contrasts within each series ============================
nboot <- as.integer(opt("n_bootstrap", 200))
min_change <- as.numeric(opt("min_log2_change", 0.5))
min_sites <- 5L
T12 <- NULL
if (length(sel) && !is.null(sv) && sv %in% names(S) && !is.null(tv) && tv %in% names(S) && nboot > 0) {
  ser <- S[!is.na(get(sv)) & get(sv) != "" & !is.na(get(tv)), .(sample_id, series = as.character(get(sv)), time = get(tv))]
  ser <- ser[series %in% ser[, uniqueN(time), by = series][V1 >= 2, series]]
  if (nrow(ser)) {
    boot_one <- function(sid) {       # log2 (peak density / background density) per family, with bootstrap draws over windows
      w <- W[sample_id == sid, .(window_id, peak, width)]
      s <- dcast(SI[sample_id == sid], window_id ~ family_id, value.var = "n_sites", fill = 0)
      M <- matrix(0, nrow(w), length(sel), dimnames = list(NULL, sel))
      if (nrow(s)) { i <- match(s$window_id, w$window_id); ok <- !is.na(i); cn <- intersect(sel, names(s)); M[i[ok], cn] <- as.matrix(s[ok, ..cn]) }
      dens <- function(wt, rows) (as.numeric(crossprod(wt, M[rows, , drop = FALSE])) + 0.5) / (sum(wt * w$width[rows]) / 1000)
      pr <- which(w$peak == 1); br <- which(w$peak == 0)
      if (!length(pr) || !length(br)) return(NULL)
      obs <- log2(dens(rep(1, length(pr)), pr) / dens(rep(1, length(br)), br))
      bs <- vapply(seq_len(nboot), function(i) log2(dens(tabulate(sample.int(length(pr), replace = TRUE), length(pr)), pr) /
                                                   dens(tabulate(sample.int(length(br), replace = TRUE), length(br)), br)), numeric(length(sel)))
      list(obs = setNames(obs, sel), bs = matrix(bs, nrow = length(sel), dimnames = list(sel, NULL)),
           n_pk = setNames(colSums(M[pr, , drop = FALSE]), sel), n_bg = setNames(colSums(M[br, , drop = FALSE]), sel))
    }
    B <- parallel::mclapply(ser$sample_id, boot_one, mc.cores = threads); names(B) <- ser$sample_id
    T12 <- rbindlist(lapply(unique(ser$series), function(g) {
      x <- ser[series == g][order(time)]; ref <- x$sample_id[1]
      if (is.null(B[[ref]])) return(NULL)
      rbindlist(lapply(seq_len(nrow(x)), function(i) { b <- B[[x$sample_id[i]]]; if (is.null(b)) return(NULL)
        d <- b$bs - B[[ref]]$bs; r0 <- B[[ref]]
        se <- apply(d, 1, sd)
        ## a contrast is tested only if both libraries have motif sites (otherwise both densities are the pseudocount
        ## and the "difference" only reflects how much sequence each library has) and the bootstrap varies
        testable <- i > 1 & (b$n_pk + b$n_bg) >= min_sites & (r0$n_pk + r0$n_bg) >= min_sites & is.finite(se) & se > 0
        data.table(series = g, family_id = sel, time = x$time[i], sample_id = x$sample_id[i], reference_time = x$time[1],
                   n_sites_peaks = b$n_pk, n_sites_background = b$n_bg,
                   log2_density_ratio = b$obs, ci_low = apply(b$bs, 1, quantile, 0.025), ci_high = apply(b$bs, 1, quantile, 0.975),
                   contrast_vs_reference = if (i == 1) NA_real_ else b$obs - r0$obs,
                   contrast_ci_low = if (i == 1) NA_real_ else apply(d, 1, quantile, 0.025), contrast_ci_high = if (i == 1) NA_real_ else apply(d, 1, quantile, 0.975),
                   contrast_se = if (i == 1) NA_real_ else se,
                   p_value = ifelse(testable, 2 * pnorm(-abs(b$obs - r0$obs) / se), NA_real_)) }))
    }))
    if (!is.null(T12) && nrow(T12)) {
      T12[, q_value := bh(p_value), by = series]
      T12[, differs_from_reference := !is.na(q_value) & q_value < alpha & abs(contrast_vs_reference) >= min_change]
      T12[, pwm_consensus := FS$pwm_consensus[match(family_id, FS$family_id)]]
      setcolorder(T12, c("series", "family_id", "pwm_consensus"))
      wt(T12, cfg, "table12_timepoint_contrasts")
      tlab <- cov_lab(tv)
      for (g in unique(T12$series)) {
        x <- T12[series == g]
        rk <- x[!is.na(p_value), .(d = any(differs_from_reference), q = min(q_value), eff = max(abs(contrast_vs_reference))), by = family_id][order(-d, q, -eff)]
        if (!nrow(rk)) next
        nsig <- rk[d == TRUE, .N]; top <- head(rk$family_id, 20)
        x <- x[family_id %in% top]
        x[, fam := factor(fam_label(family_id), levels = fam_label(top))]
        x[, tf := factor(time, levels = sort(unique(time)))]
        p20 <- ggplot(x, aes(tf, log2_density_ratio, group = 1)) +
          geom_line(colour = "#2a78d6") + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, colour = "#2a78d6") +
          geom_point(aes(shape = factor(differs_from_reference, levels = c(FALSE, TRUE))), colour = "#2a78d6", size = 2.2, show.legend = TRUE) +
          scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), limits = c("FALSE", "TRUE"), labels = c(`FALSE` = "first time point of the series, or not different from it",
                                                                                         `TRUE` = sprintf("differs from the first time point (q < %g and change >= %g log2 units)", alpha, min_change)), name = NULL, drop = FALSE) +
          facet_wrap(~ fam, ncol = 4, scales = "free_y") +
          labs(x = paste0(tlab, " (each time point is one library)"), y = "log2 motif density, peaks / background", title = paste0("Series ", g)) + theme_x(9)
        save_plot(p20, cfg, paste0("fig20_timepoint_contrasts_", gsub("[^A-Za-z0-9.-]+", "_", g)), 13, 1.8 + 2.3 * ceiling(length(top) / 4))
        fig_note(cfg, paste0("fig20_timepoint_contrasts_", gsub("[^A-Za-z0-9.-]+", "_", g)), sprintf(
          "%d of %d testable families shown (families with a difference first, then smallest q-value, then largest change). %d families differ from the first time point at some time (q < %g and |change| >= %g log2 units). Families with fewer than %d motif sites in a library are not tested. Each time point is a single library: the intervals reflect sampling of windows within that library, not biological replication.",
          length(top), nrow(rk), nsig, alpha, min_change, min_sites))
      }
    }
  }
}
if (is.null(T12) || !nrow(T12)) summary_note("fig20", "No time-course series with at least two time points and selected motif families; no time-point contrasts were computed.")
drop_empty_notes()
msg("regression models done")
