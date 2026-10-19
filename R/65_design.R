## Stage 0 of the binding models: the sample-level design, checked before any model is fitted.
## Samples: every sample with an IP BAM and a control BAM (binding needs both); the others are listed and left out.
## Covariates: categorical_covariates plus time_var binned at binding_models.time_bins. The run / batch variable is never a
## model term; it appears only in the association table, so confounding with run is visible.
## Outputs
##   tables/design_samples.tsv          every manifest sample: included in the models or not, and why
##   tables/design_estimability.tsv     one row per covariate level: reference, estimable, or dropped (aliased with ...)
##   tables/design_associations.tsv     pairwise Cramer's V between covariates (and run)
##   work/models/estimable_terms.tsv    term, covariate, level, reference: the columns the binding models will use
##   work/models/design_matrix.tsv      sample_id + the estimable 0/1 columns
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "utils.R"))
cfg <- args_config()
bm <- cfg$binding_models %||% list()
if (isFALSE(bm$enabled)) { msg("binding models disabled (binding_models.enabled: false)"); quit(status = 0) }
MD <- file.path(cfg$work, "models"); dir.create(MD, showWarnings = FALSE, recursive = TRUE)

S0 <- samples_dt(cfg)
has <- function(v) if (v %in% names(S0)) !is.na(S0[[v]]) & nzchar(S0[[v]]) else rep(FALSE, nrow(S0))
S0[, included := has("bam_file") & has("control_bam")]
S0[, reason := fifelse(included, "", fifelse(!has("bam_file"), "no IP BAM", "no control BAM"))]
wt(S0[, .(sample_id, included, reason)], cfg, "design_samples")
S <- S0[included == TRUE]
msg(sprintf("design: %d manifest samples, %d with IP and control BAMs, %d left out%s", nrow(S0), nrow(S), nrow(S0) - nrow(S),
            if (nrow(S0) > nrow(S)) paste0(" (", paste(S0[included == FALSE, paste0(sample_id, ": ", reason)], collapse = "; "), ")") else ""))
if (nrow(S) < 2) stop("fewer than two samples have both an IP and a control BAM; the binding models need at least two")

bv <- cfg$batch_var; if (is.null(bv) || !nzchar(bv) || !bv %in% names(S)) bv <- NULL
tv <- cfg$time_var
cat_v <- intersect(unlist(cfg$categorical_covariates), names(S))
cov_lab <- function(v) {
  if (v == "time_bin") return(paste0(cov_lab(tv), " (binned)"))
  if (!is.null(bv) && v == bv) return("Pipeline run (not a model term)")
  l <- cfg$covariate_labels[[v]]; if (!is.null(l)) return(l)
  gsub("_", " ", v)
}

## ordered time categories shared across series
tb <- sort(as.numeric(unlist(bm$time_bins %||% c(2, 6, 12))))
if (!is.null(tv) && tv %in% names(S) && any(!is.na(S[[tv]]))) {
  labs <- c(paste0("<=", tb[1]), if (length(tb) > 1) paste0(">", head(tb, -1), "-", tb[-1]), paste0(">", tail(tb, 1)))
  S[, time_bin := as.character(cut(get(tv), c(-Inf, tb, Inf), labels = labs))]
  S[is.na(time_bin), time_bin := "no time"]
  S[, time_bin := factor(time_bin, levels = c(labs, "no time")[c(labs, "no time") %in% time_bin])]
}
for (v in cat_v) S[, (v) := factor(fifelse(is.na(get(v)) | get(v) == "", "NA", as.character(get(v))))]
covs <- c(cat_v, if ("time_bin" %in% names(S)) "time_bin")
single <- covs[vapply(covs, function(v) nlevels(droplevels(S[[v]])) < 2, logical(1))]
if (length(single)) msg("covariates with one level among the modelled samples (no effect estimable): ", paste(single, collapse = ", "))
covs <- setdiff(covs, single)
if (!is.null(bv)) S[, (bv) := factor(get(bv))]

## treatment-contrast design; a column is kept only if it adds rank to the intercept and the earlier columns
mm_of <- function(d, v) { m <- model.matrix(~ x, data.frame(x = droplevels(d[[v]])))[, -1, drop = FALSE]
                          colnames(m) <- paste0(v, "__", levels(droplevels(d[[v]]))[-1]); m }
build_design <- function(d, vars) {
  cols <- list(); owner <- character()
  for (v in vars) { m <- mm_of(d, v); cols[[v]] <- m; owner <- c(owner, rep(v, ncol(m))) }
  if (!length(cols)) return(list(Z = matrix(0, nrow(d), 0), owner = character(), dropped = data.table()))
  Z <- do.call(cbind, cols); keep <- rep(TRUE, ncol(Z)); dropped <- list()
  for (j in seq_len(ncol(Z))) {
    A <- cbind(1, Z[, which(keep & seq_len(ncol(Z)) < j), drop = FALSE])
    if (qr(cbind(A, Z[, j]))$rank == qr(A)$rank) {
      keep[j] <- FALSE
      co <- qr.coef(qr(A), Z[, j]); names(co)[1] <- "intercept"; co <- co[!is.na(co) & abs(co) > 1e-8]
      dropped[[length(dropped) + 1]] <- data.table(term = colnames(Z)[j], covariate = owner[j],
        aliased_with = if (length(co)) sub("^\\+", "", paste(sprintf("%+g x %s", signif(co, 3), names(co)), collapse = " ")) else "zero (constant)")
    }
  }
  list(Z = Z[, keep, drop = FALSE], owner = owner[keep], dropped = rbindlist(dropped))
}
DES <- build_design(S, covs); Z <- DES$Z
rk <- qr(cbind(1, Z))$rank
if (rk != ncol(Z) + 1) stop("internal: the kept design is not full rank")   # cannot happen with the selection above
resid_df <- nrow(S) - rk

## one row per covariate level
est <- rbindlist(lapply(covs, function(v) {
  lv <- levels(droplevels(S[[v]]))
  d <- data.table(covariate = v, covariate_label = cov_lab(v), level = lv,
                  n_samples = as.integer(table(droplevels(S[[v]]))[lv]), term = c(NA, paste0(v, "__", lv[-1])))
  d[, status := fifelse(is.na(term), "reference level", fifelse(term %in% colnames(Z), "estimable", "dropped (aliased)"))]
  d[, aliased_with := if (nrow(DES$dropped)) DES$dropped$aliased_with[match(term, DES$dropped$term)] else NA_character_]
  d[]
}))
if (length(single)) est <- rbind(est, rbindlist(lapply(single, function(v) data.table(
  covariate = v, covariate_label = cov_lab(v), level = levels(droplevels(S[[v]])), n_samples = nrow(S), term = NA_character_,
  status = "single level (not estimable)", aliased_with = NA_character_))), fill = TRUE)
wt(est[, !"term"], cfg, "design_estimability")

## pairwise associations, run included so its confounding with the covariates is visible
cramers <- function(a, b) { t <- table(droplevels(a), droplevels(b)); k <- min(dim(t)) - 1
  if (k < 1) return(NA_real_); sqrt(unname(suppressWarnings(chisq.test(t, correct = FALSE)$statistic)) / (sum(t) * k)) }
av <- c(covs, bv)
assoc <- rbindlist(if (length(av) >= 2) combn(av, 2, simplify = FALSE, FUN = function(p) {
  v <- cramers(S[[p[1]]], S[[p[2]]])
  data.table(variable_a = p[1], variable_b = p[2], label_a = cov_lab(p[1]), label_b = cov_lab(p[2]), cramers_v = signif(v, 3),
             n_samples = nrow(S),
             note = fifelse(!is.na(v) & v >= 0.99, "aliased: effects cannot be separated",
                    fifelse(!is.na(v) & v >= 0.7, "strongly associated: interpret effects with caution", "")))
}) else list())
if (!nrow(assoc)) assoc <- data.table(variable_a = character(), variable_b = character(), label_a = character(), label_b = character(),
                                      cramers_v = numeric(), n_samples = integer(), note = character())
wt(assoc[order(-cramers_v)], cfg, "design_associations")

## for the binding models
terms <- est[status == "estimable", .(term, covariate, level, reference = vapply(covariate, function(v) levels(droplevels(S[[v]]))[1], ""))]
fwrite(terms, file.path(MD, "estimable_terms.tsv"), sep = "\t")
fwrite(cbind(data.table(sample_id = S$sample_id), as.data.table(Z)), file.path(MD, "design_matrix.tsv"), sep = "\t")
if (nrow(DES$dropped)) msg("aliased levels dropped: ", paste(sprintf("%s (= %s)", DES$dropped$term, DES$dropped$aliased_with), collapse = "; "))
msg(sprintf("design: %d samples, %d covariates (%s), %d estimable terms, rank %d, %d residual df",
            nrow(S), length(covs), paste(covs, collapse = ", "), ncol(Z), rk, resid_df))
if (resid_df < 1) msg("WARNING: no residual degrees of freedom at the sample level; the meta-regression cannot estimate between-sample variance")
